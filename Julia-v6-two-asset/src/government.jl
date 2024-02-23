
function tax_payment_in_dkk(LI, PC)

    # LI: Labor income
    # PC: Private contribution to retirement account

    D1 = 0.1065             # D1: Deduction 1 (deduction for people with a job): 10.65% 
    CAP = 43500.0           # CAP: max deduction 43,500DKK
    LMC = 0.08              # LMC: Labor market contribution: 8%
    D2 = 47400.0            # D2: Deduction 2 (called the ‘Basic deduction’): 47,400DKK
    MT = 0.256              # MT: municipal tax: 25.6%
    BT = 0.1209             # BT: Bottom tax: 12,09%
    TT = 0.15               # TT: Toptax: 15%
    TT_Threshold = 552000.0 # TT Threshold: 552,000 DKK
    # PC_Threshold            # I think this is the cap on the deductability of pension contributions, but not sure -- to check with others

    # Yes, contributions to annuity and life annuity accounts are tax deductible. For employer accounts it is deducted at the payroll. For private accounts, the bank/pension company reports it to the tax authorities such that you get your deduction at when the yearly summary is completed (typically March in the year following the tax year). There is also a type of pension account called Alderspension where there is no deduction, but it is then not taxed when paid out.

    # Fixed term annuities have a cap at 60900 (2023) and the cap covers the sum of contributions to employer and private accounts. Life annuities in employer accounts have no cap, but for private life annuity accounts there is a cap at 56.100. The cap on Alderspension is 8800. (There is a higher cap at 56900 if you are within 7 years from your social security eligibility age.)

    # So, for employer accounts it is already deducted from the amount of earnings that we see in the admin data, i.e. LI is already net of employer contributions. For private contributions it would be deductible . For private contributions you first pay Labour Market Contribution (LMC) but contributions are then deductible in the municipal tax (MT) and the bottom (BT) and the top tax (TT). I have tried to insert this as P in the equation below. It does not consider all the details but is perhaps a reasonable way to approximate the system.

    # tax_payment = LI*(1.0 - max(D1*LI, CAP)) * LMC + 
    #               LI*(1.0 - LMC - D2 - PC*1[PC<PC_Threshold])* (MT+BT) + 
    #               LI*(1.0 - LMC - D2 - PC*1[PC<PC_Threshold]))* TT*1[LI*(1 - LMC - D2)>TT_Threshold]

                  # NOTE: originally included PC*1[PC<PC_Threshold]. Should this be rewritten as min(PC, PC_Threshold) ?

    indicator_above_TT_Threshold = (LI*(1.0 - LMC - D2)>TT_Threshold)

    tax_payment = LI*(1.0 - max(D1*LI, CAP)) * LMC + 
                  LI*(1.0 - LMC - D2)* (MT+BT) + 
                  LI*(1.0 - LMC - D2)* TT * indicator_above_TT_Threshold

    return tax_payment

end

LIs = collect(0.0:100000.:700000.)

tax = [ tax_payment_in_dkk(LI, 0.0) for LI in LIs]

plot(LIs./1000, tax./1000)