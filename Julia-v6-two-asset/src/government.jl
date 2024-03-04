using Plots

function tax_payment_in_dkk(LI, PC)

    # LI: Labor income
    # PC: Private pension contribution

    # https://themis.dk/synopsis/docs/Oversigt_over_beloebsgraenser.html
    D1 = 0.044              # D1: Deduction 1 (Beskæftigelsesfradrag, in percent)
    D1max = 14100.0         # D1max: max deduction (Beskæftigelsesfradrag, max)
    LMC = 0.08              # LMC: Labor market contribution
    D2 = 42900.0            # D2: Deduction 2 (Personfradrag)
    MT = 0.253              # MT: municipal tax
    BT = 0.1164             # BT: Bundskat + Sundhedsbidrag
    TT = 0.15               # TT: Topskat
    TT_thld = 389900.0      # TT Threshold
    DPmax = 46000.0         # max. pension contribution

    # Yes, contributions to annuity and life annuity accounts are tax deductible. For employer accounts it is deducted at the payroll. For private accounts, the bank/pension company reports it to the tax authorities such that you get your deduction at when the yearly summary is completed (typically March in the year following the tax year). There is also a type of pension account called Alderspension where there is no deduction, but it is then not taxed when paid out.

    # Fixed term annuities have a cap at 60900 (2023) and the cap covers the sum of contributions to employer and private accounts. Life annuities in employer accounts have no cap, but for private life annuity accounts there is a cap at 56.100. The cap on Alderspension is 8800. (There is a higher cap at 56900 if you are within 7 years from your social security eligibility age.)

    # So, for employer accounts it is already deducted from the amount of earnings that we see in the admin data, i.e. LI is already net of employer contributions. For private contributions it would be deductible . For private contributions you first pay Labour Market Contribution (LMC) but contributions are then deductible in the municipal tax (MT) and the bottom (BT) and the top tax (TT). I have tried to insert this as P in the equation below. It does not consider all the details but is perhaps a reasonable way to approximate the system.

    DP = min(PC, DPmax)
    tax_base_lmc = LI - DP - min(D1*(LI-DP), D1max)
    ind_top = tax_base_lmc * (1-LMC) - D2 > TT_thld
    
    tax_payment =   tax_base_lmc * LMC + 
                    max(tax_base_lmc * (1-LMC) - D2, 0) * (MT + BT) + 
                    (tax_base_lmc * (1-LMC) - D2 - TT_thld) * TT * ind_top

    return tax_payment


end

LIs = collect(0.0:10000.:700000.)

tax = [tax_payment_in_dkk(LI, 0.0) for LI in LIs]
ATR = [tax_payment_in_dkk(LI, 0.0)/LI for LI in LIs]

plot(LIs./1000, tax./1000, title="taxes", legend = false)
plot(LIs./1000, ATR./1000, title="average tax rate", legend = false)