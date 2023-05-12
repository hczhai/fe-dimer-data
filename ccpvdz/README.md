
## HF / DFT

* mf-1        UHF
* mf-2        UKS-TPSS
* mf-3        UKS-B3LYP

## CCSD / CCSD(T)

* cc-5        UHF/CCSD(T)
* cc-6        UKS-TPSS/CCSD(T)
* cc-7        UKS-B3LYP/CCSD(T)
* cc-13       UHF/CCSD <S^2>

## Active space (based on UHF/CCSD natural orbitals)

* select 17   (36o, 48e)
* select 18   (55o, 48e)
* select 19   (63o, 64e)
* select 20   (88o, 64e)

## Active space DMRG / CCSD / CCSD(T)

* casci 21/22/23/24     active space DMRG
* casci 25/26/27/28     active space CCSD
* casci 29/30/31/32     active space CCSD(T)

## FCIDUMP files

* for: casci 21/22/23/24
* decompression: ``tar -Jxvf FCIDUMP-24.xz``
