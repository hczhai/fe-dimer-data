
"""Parsing outputs to generate main.xlsx."""

import os

configs = ["HC", "HS", "HFe", "HFe2"]

def get_data(dir, runs):

    data = {}
    pwd = os.popen('pwd').read().strip()
    os.chdir("../" + dir)

    for config in configs:
        data[config] = {}
        os.chdir('%s' % config.lower())
        for run in runs:
            opt = os.popen('hife log %s exact' % run).read().strip()
            try:
                if "mf " in run:
                    if "NO CONV" in opt:
                        ener = opt.split("\n")[-1].split("NITER =")[0].split("NO CONV !!!")[1].strip()
                        ssq = "0.0"
                    else:
                        ener = opt.split("\n")[-1].split("NITER =")[0].split("E =")[1].strip()
                        ssq = opt.split("\n")[0].split("S^2 =")[-1].strip()
                    f = opt.split("\n")[-2].split("FILE =")[1].split("JOBID")[0].strip()
                    fespin = [0.0, 0.0]
                    with open(f, "r") as f:
                        for l in f.readlines():
                            if l.startswith("charge of   20Fe"):
                                fespin[0] = float(l.split()[-1].strip())
                            elif l.startswith("charge of   21Fe"):
                                fespin[1] = float(l.split()[-1].strip())
                    data[config][run] = [float(ener), float(ssq), "(%.1f, %.1f)" % tuple(fespin)]
                elif "cc " in run:
                    ener = opt.split("\n")[-2].split("ECCSD =")[1].split("ECCSD(T)")[0].strip()
                    ener_t = opt.split("\n")[-2].split("ECCSD(T) =")[1].strip()
                    data[config][run] = [float(ener), float(ener_t)]
                elif "casci " in run and "DMRG" in opt:
                    dmrg_fn = opt.split("\n")[-2].split("FILE =")[1].split()[0]
                    with open(dmrg_fn, 'r') as df:
                        for dfl in df.readlines():
                            if "EXTRAP Energy" in dfl:
                                ener, ener_err = dfl.split()[3], dfl.split()[5]
                                break
                    csf_fn = dmrg_fn.split("-rev.out")[0] + "-csf.out.1"
                    with open(csf_fn, 'r') as df:
                        for dfl in df.readlines():
                            if "CSF          0" in dfl:
                                csf = dfl.split()[4]
                                break
                    data[config][run] = [float(ener), float(ener_err), float(csf) ** 2]
                elif "casci " in run:
                    ener = opt.split("\n")[-1].split("CASCI E =")[1].split()[0]
                    data[config][run] = [float(ener)]
                elif "mp " in run:
                    ener = opt.split("\n")[-1].split("EMP2 =")[1].split("NITER")[0].strip()
                    data[config][run] = [float(ener)]
            except IndexError as x:
                print("I CONFIG = %s RUN = %s LOG = " % (config, run))
                print(opt)
                raise x
            except ValueError as x:
                print("V CONFIG = %s RUN = %s LOG = " % (config, run))
                print(opt)
                raise x
        os.chdir("..")

    os.chdir(pwd)
    return data

def write_data(worksheet, runs, prop, start_row, start_col, e_format, data, data_idx,
    prev_row=0, write_header=True, write_left=True, is_num=True, col_shift=0, unit_scale=False):

    if write_left:
        worksheet.write(start_row + write_header, start_col, prop)
        start_col += 1

    from xlsxwriter.utility import xl_rowcol_to_cell as XRC

    if write_left:
        if data_idx != "average":
            for ip, p in enumerate(configs):
                worksheet.write(start_row + ip + write_header, start_col, p)
        else:
            worksheet.write(start_row + 0 + write_header, start_col, "average")
        start_col += 1

    if write_header:
        for ir, r in enumerate(runs):
            worksheet.write(start_row, start_col + ir, runs_map[r])
        start_row += 1

    start_col += col_shift
    
    if data_idx == "delta_e":
        for ir, r in enumerate(runs):
            for ip, p in enumerate(configs):
                worksheet.write_formula(start_row + ip, start_col + ir,
                    "=(%s-%s)*%s" % (
                        XRC(prev_row + 1 + ip, start_col + ir),
                        XRC(prev_row + 1 + 0, start_col + ir, True, False),
                        XRC(0, 1, True, True),
                    ), e_format)
    elif data_idx == "de_average":
        for ir, r in enumerate(runs):
            for ip, p in enumerate(configs):
                worksheet.write_formula(start_row + ip, start_col + ir,
                    "=(%s-AVERAGE(%s:%s))*%s" % (
                        XRC(prev_row + 1 + ip, start_col + ir),
                        XRC(prev_row + 1 + 0, start_col + ir, True, False),
                        XRC(prev_row + 1 + 3, start_col + ir, True, False),
                        XRC(0, 1, True, True),
                    ), e_format)
    elif data_idx == "average":
        for ir, r in enumerate(runs):
            worksheet.write_formula(start_row + 0, start_col + ir,
                    "=AVERAGE(%s:%s)" % (
                        XRC(prev_row + 1 + 0, start_col + ir, True, False),
                        XRC(prev_row + 1 + 3, start_col + ir, True, False),
                    ), e_format)
    elif data_idx == "e_corr":
        for ir, r in enumerate(runs):
            for ip, p in enumerate(configs):
                worksheet.write_formula(start_row + ip, start_col + ir,
                    "=(%s-%s)*%s" % (
                        XRC(prev_row + 1 + ip, start_col + ir),
                        XRC(prev_row + 1 + ip, start_col + 0, True, False),
                        XRC(0, 1, True, True),
                    ), e_format)
    elif unit_scale and is_num:
        for ir, r in enumerate(runs):
            for ip, p in enumerate(configs):
                worksheet.write_formula(start_row + ip, start_col + ir, "=%20.15f*%s" % (
                    data[p][r][data_idx], XRC(0, 1, True, True)), e_format)
    else:
        for ir, r in enumerate(runs):
            for ip, p in enumerate(configs):
                if not is_num:
                    worksheet.write(start_row + ip, start_col + ir, data[p][r][data_idx], e_format)
                else:
                    worksheet.write_number(start_row + ip, start_col + ir, data[p][r][data_idx], e_format)

import xlsxwriter

workbook = xlsxwriter.Workbook('main.xlsx')
e_format = workbook.add_format({'num_format': '###0.000000'})
eav_format = workbook.add_format({'num_format': '###0.0000'})
de_format = workbook.add_format({'num_format': '###0.0'})
dde_format = workbook.add_format({'num_format': '###0.00'})
s_format = workbook.add_format()
s_format.set_align('right')

### DFT functionals

worksheet = workbook.add_worksheet('SCF')
worksheet.write(0, 0, 'Hartree to kJ/mol')
worksheet.write(0, 1, 2625.5)
worksheet.set_column(0, 0, len('Hartree to kJ/mol'))

start_col = 0
write_left = True

runs = ["mf 6", "mf 10", "mf 11", "mf 23", "mf 39", "mf 21", "mf 35", "mf 7", "mf 33", "mf 22"]
runs_map = {
    "mf 6": "UKS/TPSS", "mf 10": "UKS/BLYP", "mf 11": "UKS/PBE",
    "mf 23": "UKS/B97-D", "mf 39": "UKS/r2SCAN", "mf 21": "UKS/TPSSh",
    "mf 35": "UKS/B3LYP*", "mf 7": "UKS/B3LYP", "mf 33": "UKS/PBE0",
    "mf 22": "UKS/M06",
}

start_row = 2
data = get_data("ccpvqz", runs)
worksheet.write(start_row, start_col + 2, 'ccpvqz-dk-x2c')

start_row += 2
write_data(worksheet, runs, 'Energy in Hartree', start_row, start_col, e_format,
    data, data_idx=0, write_left=write_left)

prev_row = start_row
start_row += len(configs)  + 2

write_data(worksheet, runs, 'Delta E in kJ/mol', start_row, start_col, de_format, data, data_idx='delta_e',
    prev_row=prev_row, write_header=False, write_left=write_left)
start_row += len(configs)  + 1

worksheet.set_column(2, worksheet.dim_colmax, len('######0.000000'))

### CC Energies in DZ

worksheet = workbook.add_worksheet('CC')
worksheet.write(0, 0, 'Hartree to kJ/mol')
worksheet.write(0, 1, 2625.5)
worksheet.set_column(0, 0, len('Hartree to kJ/mol'))

start_col = 0
write_left = True

runs = ["mf 1", "mf 2", "mf 3", "cc 5", "cc 6", "cc 7"]
runs_map = {
    "mf 1": "UHF", "mf 2": "UKS-TPSS", "mf 3": "UKS-B3LYP",
    "cc 5": "UHF/CCSD", "cc 6": "UKS-TPSS/CCSD", "cc 7": "UKS-TPSS/CCSD",
    "cc 5*": "UHF/CCSD(T)", "cc 6*": "UKS-TPSS/CCSD(T)", "cc 7*": "UKS-TPSS/CCSD(T)",
}

start_row = 2
data = get_data("ccpvdz", runs)
for k, v in data.items():
    for kk in list(v):
        if 'cc' in kk:
            v[kk + "*"] = v[kk][1:]
runs += ["cc 5*", "cc 6*", "cc 7*"]
worksheet.write(start_row, start_col + 2, 'ccpvdz-dk-x2c')

start_row += 2
write_data(worksheet, runs, 'Energy in Hartree', start_row, start_col, e_format,
    data, data_idx=0, write_left=write_left)

prev_row = start_row
start_row += len(configs)  + 2

write_data(worksheet, runs, 'Delta E in kJ/mol', start_row, start_col, de_format, data, data_idx='delta_e',
    prev_row=prev_row, write_header=False, write_left=write_left)
start_row += len(configs)  + 1

worksheet.set_column(2, worksheet.dim_colmax, len('######0.000000'))

### Active space in DZ

worksheet = workbook.add_worksheet('active')
worksheet.write(0, 0, 'Hartree to kJ/mol')
worksheet.write(0, 1, 2625.5)
worksheet.set_column(0, 0, len('Hartree to kJ/mol'))

start_col = 0
write_left = True

col_names = ["(36o, 48e)", "(55o, 48e)", "(63o, 64e)", "(88o, 64e)", "full space"]

for ii in range(len(col_names)):

    start_row = 2

    if col_names[ii] == "full space":
        runs_active = ["cc 5"]
        runs_map = { "cc 5": "CASCI/CCSD", "cc 5*": "CASCI/CCSD(T)" }
        data = get_data("ccpvdz", runs_active)
        for k, v in data.items():
            for kk in list(v):
                if 'cc' in kk:
                    v[kk + "*"] = v[kk][1:]
        runs_active += ["cc 5*"]
    else:
        runs_active = ["casci %d" % (25 + ii), "casci %d" % (29 + ii), "casci %d" % (21 + ii)]
        runs_map = {
            "casci %d" % (25 + ii): "CASCI/CCSD",
            "casci %d" % (29 + ii): "CASCI/CCSD(T)",
            "casci %d" % (21 + ii): "CASCI/DMRG",
        }
        data = get_data("ccpvdz", runs_active)

    worksheet.write(start_row, start_col + (2 if ii == 0 else 0), col_names[ii])

    start_row += 2
    write_data(worksheet, runs_active, 'Energy in Hartree', start_row, start_col, e_format,
        data, data_idx=0, write_left=write_left)

    prev_row = start_row
    start_row += len(configs)  + 2

    write_data(worksheet, runs_active, 'Delta E in kJ/mol', start_row, start_col, de_format, data, data_idx='delta_e',
        prev_row=prev_row, write_header=False, write_left=write_left)

    start_row += len(configs)  + 2

    write_data(worksheet, runs_active, 'E - E[CCSD] in kJ/mol', start_row, start_col, de_format, data, data_idx='e_corr',
        prev_row=prev_row, write_header=False, write_left=write_left)

    if col_names[ii] != "full space":
        start_row += len(configs)  + 2

        write_data(worksheet, runs_active[-1:], 'DMRG extra. error in kJ/mol', start_row, start_col, de_format, data, data_idx=1,
            write_header=False, write_left=write_left, col_shift=2, unit_scale=True)

        start_row += len(configs)  + 2

        write_data(worksheet, runs_active[-1:], 'Max CSF weight', start_row, start_col, dde_format, data, data_idx=2,
            write_header=False, write_left=write_left, col_shift=2)

    start_row += len(configs)  + 1

    start_col += 3 if ii != 0 else 5
    write_left = False

worksheet.set_column(2, worksheet.dim_colmax, len('######0.000000'))

### Basis sets

worksheet = workbook.add_worksheet('Basis set')
worksheet.write(0, 0, 'Hartree to kJ/mol')
worksheet.write(0, 1, 2625.5)
worksheet.set_column(0, 0, len('Hartree to kJ/mol'))

data_svp = get_data("def2-svp", ["mf 1", "cc 5", "mp 17"])
data_dz = get_data("ccpvdz", ["mf 1", "cc 5"])
data_tz = get_data("ccpvtz", ["mf 1", "cc 8"])
data_qz = get_data("ccpvqz", ["mf 1"])
data_mp2 = get_data("mp2", ["mp 17", "mp 21", "mp 25"])

runs_mf = ["mf 1s", "mf 1d", "mf 1t", "mf 1q"]
runs_mp2 = ["mp 17s", "mp 17d", "mp 21t", "mp 25q", "mp tql", "mp dtl"]
runs_cc = ["cc 5s", "cc 5d", "cc 8t", "cc dtl"]
runs_cc_t = ["cct 5s", "cct 5d", "cct 8t", "cct dtl"]
data = { k : {} for k in data_svp }
for dd, ds in zip([data_svp, data_dz, data_tz, data_qz], "sdtq"):
    for k, v in dd.items():
        for kk, vv in v.items():
            data[k][kk + ds] = vv
for k, v in data_mp2.items():
    for kk, vv in v.items():
        if kk == "mp 17":
            data[k][kk + 'd'] = vv
        elif kk == "mp 21":
            data[k][kk + 't'] = vv
        elif kk == "mp 25":
            data[k][kk + 'q'] = vv

for k, v in data.items():
    for kk, vv in list(v.items()):
        if "cc" in kk:
            v["cct " + kk.split()[1]] = [vv[1] - vv[0]]
    for kk, vv in list(v.items()):
        if "mf" not in kk and "cct" not in kk and kk != "mp 25q":
            vv[0] = vv[0] - v["mf 1" + kk[-1]][0]
    mp2_d_corr = data[k]["mp 17d"][0]
    mp2_t_corr = data[k]["mp 21t"][0]
    mp2_q_corr = data[k]["mp 25q"][0]
    cc_d_corr = data[k]["cc 5d"][0]
    cc_t_corr = data[k]["cc 8t"][0]
    cct_d_corr = data[k]["cct 5d"][0]
    cct_t_corr = data[k]["cct 8t"][0]
    data[k]["mp dtl"] = [(2 ** 2.4 * mp2_d_corr - 3 ** 2.4 * mp2_t_corr) / (2 ** 2.4 - 3 ** 2.4)]
    data[k]["mp tql"] = [(3 ** 2.4 * mp2_t_corr - 4 ** 2.4 * mp2_q_corr) / (3 ** 2.4 - 4 ** 2.4)]
    data[k]["cc dtl"] = [(2 ** 2.4 * cc_d_corr - 3 ** 2.4 * cc_t_corr) / (2 ** 2.4 - 3 ** 2.4)]
    data[k]["cct dtl"] = [(2 ** 2.4 * cct_d_corr - 3 ** 2.4 * cct_t_corr) / (2 ** 2.4 - 3 ** 2.4)]

start_col = 0
write_left = True

col_names = ["UHF", "UHF/MP2 (corr)", "UHF/CCSD (corr)", "UHF/CCSD(T) [(T) only]"]

for ii in range(len(col_names)):

    start_row = 2

    if col_names[ii] == "UHF":
        runs = runs_mf
        runs_map = { "mf 1s": "def2-SV(P)", "mf 1d": "cc-pVDZ-DK", "mf 1t": "cc-pVTZ-DK", "mf 1q": "cc-pVQZ-DK" }
    elif col_names[ii] == "UHF/MP2 (corr)":
        runs = runs_mp2
        runs_map = { "mp 17s": "def2-SV(P)", "mp 17d": "cc-pVDZ-DK", "mp 21t": "cc-pVTZ-DK",
            "mp 25q": "cc-pVQZ-DK", "mp tql": "TZ/QZ CBS", "mp dtl": "DZ/TZ CBS" }
    elif col_names[ii] == "UHF/CCSD (corr)":
        runs = runs_cc
        runs_map = { "cc 5s": "def2-SV(P)", "cc 5d": "cc-pVDZ-DK", "cc 8t": "cc-pVTZ-DK",
            "cc dtl": "DZ/TZ CBS" }
    elif col_names[ii] == "UHF/CCSD(T) [(T) only]":
        runs = runs_cc_t
        runs_map = { "cct 5s": "def2-SV(P)", "cct 5d": "cc-pVDZ-DK", "cct 8t": "cc-pVTZ-DK",
            "cct dtl": "DZ/TZ CBS" }

    worksheet.write(start_row, start_col + (2 if ii == 0 else 0), col_names[ii])

    start_row += 2
    write_data(worksheet, runs, 'E[corr] in Hartree', start_row, start_col, e_format,
        data, data_idx=0, write_left=write_left)

    prev_row = start_row
    start_row += len(configs) + 2

    write_data(worksheet, runs, 'Average E[corr]', start_row, start_col, eav_format, data, data_idx='average',
        prev_row=prev_row, write_header=False, write_left=write_left)

    start_row += 1 + 2

    write_data(worksheet, runs, 'E[corr] - E[av] in kJ/mol', start_row, start_col, de_format, data, data_idx='de_average',
        prev_row=prev_row, write_header=False, write_left=write_left)

    start_row += len(configs)  + 1

    start_col += len(runs) if ii != 0 else 2 + len(runs)
    write_left = False

worksheet.set_column(2, worksheet.dim_colmax, len('######0.000000'))

### END

workbook.close()