/* /////////////////////////////////////////////////////////////////// */
   ////////////////////// erf coefficients /////////////////////////

`define  ERF_A1     3.16112374387056560e0
`define  ERF_A2     1.13864154151050156e2
`define  ERF_A3     3.77485237685302021e2
`define  ERF_A4     3.20937758913846947e3
`define  ERF_A5     1.85777706184603153e-1
`define  ERF_B1     2.36012909523441209e1
`define  ERF_B2     2.44024637934444173e2
`define  ERF_B3     1.28261652607737228e3
`define  ERF_B4     2.84423683343917062e3
   
`define  ERF_C1     5.64188496988670089e-1
`define  ERF_C2     8.88314979438837594e0
`define  ERF_C3     6.61191906371416295e1
`define  ERF_C4     2.98635138197400131e2
`define  ERF_C5     8.81952221241769090e2
`define  ERF_C6     1.71204761263407058e3
`define  ERF_C7     2.05107837782607147e3
`define  ERF_C8     1.23033935479799725e3
`define  ERF_C9     2.15311535474403846e-8
`define  ERF_D1     1.57449261107098347e1
`define  ERF_D2     1.17693950891312499e2
`define  ERF_D3     5.37181101862009858e2
`define  ERF_D4     1.62138957456669019e3
`define  ERF_D5     3.29079923573345963e3
`define  ERF_D6     4.36261909014324716e3
`define  ERF_D7     3.43936767414372164e3
`define  ERF_D8     1.23033935480374942e3

`define  ERF_P1     3.05326634961232344e-1
`define  ERF_P2     3.60344899949804439e-1
`define  ERF_P3     1.25781726111229246e-1
`define  ERF_P4     1.60837851487422766e-2
`define  ERF_P5     6.58749161529837803e-4
`define  ERF_P6     1.63153871373020978e-2
`define  ERF_Q1     2.56852019228982242e0
`define  ERF_Q2     1.87295284992346047e0
`define  ERF_Q3     5.27905102951428412e-1
`define  ERF_Q4     6.05183413124413191e-2
`define  ERF_Q5     2.33520497626869185e-3

`define  ERF_SQRPI  5.6418958354775628695e-1

/* /////////////////////////////////////////////////////////////////// */
   ////////////////////// erf functions /////////////////////////

analog function real f_erf_r1_sub;
    input x;
    real x;
    real xsq, num, den;
    begin
        xsq = x*x;
        num = (((`ERF_A5*xsq + `ERF_A1)*xsq + `ERF_A2)*xsq + `ERF_A3)*xsq;
        den = (((xsq + `ERF_B1)*xsq + `ERF_B2)*xsq + `ERF_B3)*xsq;
        f_erf_r1_sub = x*(num+`ERF_A4)/(den+`ERF_B4);
    end
endfunction

analog function real f_erf_r2_sub;
    input x;
    real x;
    real y, delf, num, den;
    begin
        y = floor(16.0*x)/16.0;
        delf = exp(-(y*y)-(x-y)*(x+y));
        num = (((((((`ERF_C9*x + `ERF_C1)*x + `ERF_C2)*x + `ERF_C3)*x + `ERF_C4)*x + `ERF_C5)*x + `ERF_C6)*x + `ERF_C7)*x;
        den = (((((((x + `ERF_D1)*x + `ERF_D2)*x + `ERF_D3)*x + `ERF_D4)*x + `ERF_D5)*x + `ERF_D6)*x + `ERF_D7)*x;
        f_erf_r2_sub = 1.0 - delf*(num+`ERF_C8)/(den+`ERF_D8);
    end
endfunction

analog function real f_erf_r3_sub;
    input x;
    real x;
    real y, z, delf, num, den;
    begin
        y = floor(16.0*x)/16.0;
        z = 1.0/(x*x);
        delf = exp(-(y*y)-(x-y)*(x+y));
        num = ((((`ERF_P6*z + `ERF_P1)*z + `ERF_P2)*z + `ERF_P3)*z + `ERF_P4)*z;
        den = ((((z + `ERF_Q1)*z + `ERF_Q2)*z + `ERF_Q3)*z + `ERF_Q4)*z;
        f_erf_r3_sub = 1.0 - delf*(`ERF_SQRPI - z*(num + `ERF_P5)/(den + `ERF_Q5)) / x;
    end
endfunction

analog function real f_erf_subroutine;
    input x;
    real x;
    begin
        if (x <= 0.46875)
            f_erf_subroutine = f_erf_r1_sub(x);
        else if (x <= 4.0)
            f_erf_subroutine = f_erf_r2_sub(x);
        else if (x <= 5.8)
            f_erf_subroutine = f_erf_r3_sub(x);
        else
            f_erf_subroutine = 1.0;
    end
endfunction

analog function real erf;
    input x;
    real x;
    begin
        if (x < 0.0)
            erf = -f_erf_subroutine(-x);
        else
            erf = f_erf_subroutine(x);
    end
endfunction

analog function real erfc;
    input x;
    real x;
    begin
        erfc = 1.0 - erf(x);
    end
endfunction


