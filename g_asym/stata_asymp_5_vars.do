preserve

keep if type == "FIT"

statsby ///
    b_ofi_bin1 = _b[ofi_bin1] ///
    b_ofi_bin2_pos = _b[ofi_bin2_pos] ///
    b_ofi_bin2_neg = _b[ofi_bin2_neg] ///
    b_ofi_bin3_pos = _b[ofi_bin3_pos] ///
    b_ofi_bin3_neg = _b[ofi_bin3_neg] ///
    N = e(N) ///
    r2 = e(r2), ///
    by(yyyymm) clear : ///
    regress ret ///
        ofi_bin1 ///
        ofi_bin2_pos ofi_bin2_neg ///
        ofi_bin3_pos ofi_bin3_neg ///
		accruals agr beta bm gp indmom intmom iss_1y iss_5y mom mom_season noa rev  ///
		effective_spread quoted_spread realized_vol size turnover dollar_volume
		
		
		


save "fm_monthly_fit_ofi_posneg.dta", replace

foreach x in ///
    ofi_bin1  ///
    ofi_bin2_pos ofi_bin2_neg ///
    ofi_bin3_pos ofi_bin3_neg {

    quietly summarize b_`x'
    local mean_`x' = r(mean)
    local se_`x'   = r(sd) / sqrt(r(N))
    local t_`x'    = `mean_`x'' / `se_`x''
}

display "Fama-MacBeth estimates"
display "Sample: type == FIT"
display "--------------------------------------------------"

foreach x in ///
    ofi_bin1 ///
    ofi_bin2_pos ofi_bin2_neg ///
    ofi_bin3_pos ofi_bin3_neg {

    display "`x': " ///
        %10.4f `mean_`x'' ///
        "   SE: " %10.4f `se_`x'' ///
        "   t: "  %10.3f `t_`x''
}

* Pairwise differences
gen d_pos_2_1 = b_ofi_bin2_pos - b_ofi_bin1_pos
gen d_pos_3_2 = b_ofi_bin3_pos - b_ofi_bin2_pos
gen d_pos_3_1 = b_ofi_bin3_pos - b_ofi_bin1_pos

gen d_neg_2_1 = b_ofi_bin2_neg - b_ofi_bin1_neg
gen d_neg_3_2 = b_ofi_bin3_neg - b_ofi_bin2_neg
gen d_neg_3_1 = b_ofi_bin3_neg - b_ofi_bin1_neg

foreach x in ///
    d_pos_2_1 d_pos_3_2 d_pos_3_1 ///
    d_neg_2_1 d_neg_3_2 d_neg_3_1 {

    quietly summarize `x'
    local mean_`x' = r(mean)
    local se_`x'   = r(sd) / sqrt(r(N))
    local t_`x'    = `mean_`x'' / `se_`x''
}

display " "
display "Pairwise differences (Fama-MacBeth)"
display "-----------------------------------"

display "POS:"
display "2-1: " %9.4f `mean_d_pos_2_1' "  SE: " %9.4f `se_d_pos_2_1' "  t: " %9.3f `t_d_pos_2_1'
display "3-2: " %9.4f `mean_d_pos_3_2' "  SE: " %9.4f `se_d_pos_3_2' "  t: " %9.3f `t_d_pos_3_2'
display "3-1: " %9.4f `mean_d_pos_3_1' "  SE: " %9.4f `se_d_pos_3_1' "  t: " %9.3f `t_d_pos_3_1'

display " "
display "NEG:"
display "2-1: " %9.4f `mean_d_neg_2_1' "  SE: " %9.4f `se_d_neg_2_1' "  t: " %9.3f `t_d_neg_2_1'
display "3-2: " %9.4f `mean_d_neg_3_2' "  SE: " %9.4f `se_d_neg_3_2' "  t: " %9.3f `t_d_neg_3_2'
display "3-1: " %9.4f `mean_d_neg_3_1' "  SE: " %9.4f `se_d_neg_3_1' "  t: " %9.3f `t_d_neg_3_1'

restore
