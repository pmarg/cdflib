
    ! CDFLIB - edited by Patrick Macnamara
    ! Converted to module and added kind parameters to improve portability.
    ! Removed cdf* routines which used old-fashioned Fortran code.
    ! Changed name of beta function to beta_func.
    module cdflib
    use, intrinsic :: iso_fortran_env, only: dp => real64, int32

    implicit none

    contains

    function algdiv ( a, b )

    !*****************************************************************************80
    !
    !! ALGDIV computes ln ( Gamma ( B ) / Gamma ( A + B ) ) when 8 <= B.
    !
    !  Discussion:
    !
    !    In this algorithm, DEL(X) is the function defined by
    !
    !      ln ( Gamma(X) ) = ( X - 0.5 ) * ln ( X ) - X + 0.5 * ln ( 2 * PI )
    !                      + DEL(X).
    !
    !  Author:
    !
    !    Armido DiDinato, Alfred Morris
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) A, B, define the arguments.
    !
    !    Output, real ( kind = dp ) ALGDIV, the value of ln(Gamma(B)/Gamma(A+B)).
    !
    implicit none

    real ( kind = dp ) :: a,b
    real ( kind = dp ) :: algdiv

    real ( kind = dp ) :: c,d,h,s11,s3,s5,s7,s9,t,u,v,w,x,x2

    real ( kind = dp ), parameter :: c0 =  0.833333333333333D-01
    real ( kind = dp ), parameter :: c1 = -0.277777777760991D-02
    real ( kind = dp ), parameter :: c2 =  0.793650666825390D-03
    real ( kind = dp ), parameter :: c3 = -0.595202931351870D-03
    real ( kind = dp ), parameter :: c4 =  0.837308034031215D-03
    real ( kind = dp ), parameter :: c5 = -0.165322962780713D-02

    if ( b < a ) then
        h = b / a
        c = 1.0D+00 / ( 1.0D+00 + h )
        x = h / ( 1.0D+00 + h )
        d = a + ( b - 0.5D+00 )
    else
        h = a / b
        c = h / ( 1.0D+00 + h )
        x = 1.0D+00 / ( 1.0D+00 + h )
        d = b + ( a - 0.5D+00 )
    end if
    !
    !  Set SN = (1 - X^N)/(1 - X).
    !
    x2 = x * x
    s3 = 1.0D+00 + ( x + x2 )
    s5 = 1.0D+00 + ( x + x2 * s3 )
    s7 = 1.0D+00 + ( x + x2 * s5 )
    s9 = 1.0D+00 + ( x + x2 * s7 )
    s11 = 1.0D+00 + ( x + x2 * s9 )
    !
    !  Set W = DEL(B) - DEL(A + B).
    !
    t = ( 1.0D+00 / b )**2
    w = (((( &
        c5 * s11  * t &
        + c4 * s9 ) * t &
        + c3 * s7 ) * t &
        + c2 * s5 ) * t &
        + c1 * s3 ) * t &
        + c0

    w = w * ( c / b )
    !
    !  Combine the results.
    !
    u = d * alnrel ( a / b )
    v = a * ( log ( b ) - 1.0D+00 )

    if ( v < u ) then
        algdiv = ( w - v ) - u
    else
        algdiv = ( w - u ) - v
    end if

    return
    end function algdiv

    function alnrel ( a )

    !*****************************************************************************80
    !
    !! ALNREL evaluates the function ln ( 1 + A ).
    !
    !  Author:
    !
    !    Armido DiDinato, Alfred Morris
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) A, the argument.
    !
    !    Output, real ( kind = dp ) ALNREL, the value of ln ( 1 + A ).
    !
    implicit none

    real ( kind = dp ) :: a
    real ( kind = dp ) :: alnrel

    real ( kind = dp ) :: t,t2,w,x

    real ( kind = dp ), parameter :: p1 = -0.129418923021993D+01
    real ( kind = dp ), parameter :: p2 =  0.405303492862024D+00
    real ( kind = dp ), parameter :: p3 = -0.178874546012214D-01
    real ( kind = dp ), parameter :: q1 = -0.162752256355323D+01
    real ( kind = dp ), parameter :: q2 =  0.747811014037616D+00
    real ( kind = dp ), parameter :: q3 = -0.845104217945565D-01

    if ( abs ( a ) <= 0.375D+00 ) then

        t = a / ( a +  2.0D+00  )
        t2 = t * t

        w = ((( p3 * t2 + p2 ) * t2 + p1 ) * t2 + 1.0D+00 ) &
            / ((( q3 * t2 + q2 ) * t2 + q1 ) * t2 + 1.0D+00 )

        alnrel =  2.0D+00  * t * w

    else

        x = 1.0D+00 + real ( a, kind = dp )
        alnrel = log ( x )

    end if

    return
    end function alnrel

    function apser ( a, b, x, eps )

    !*****************************************************************************80
    !
    !! APSER computes the incomplete beta ratio I(SUB(1-X))(B,A).
    !
    !  Discussion:
    !
    !    APSER is used only for cases where
    !
    !      A <= min ( EPS, EPS * B ),
    !      B * X <= 1, and
    !      X <= 0.5.
    !
    !  Author:
    !
    !    Armido DiDinato, Alfred Morris
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) A, B, X, the parameters of the
    !    incomplete beta ratio.
    !
    !    Input, real ( kind = dp ) EPS, a tolerance.
    !
    !    Output, real ( kind = dp ) APSER, the computed value of the
    !    incomplete beta ratio.
    !
    implicit none

    real ( kind = dp ) :: a,b,x,eps
    real ( kind = dp ) :: apser

    real ( kind = dp ) :: aj,bx,c,j,s,t,tol

    real ( kind = dp ), parameter :: g = 0.577215664901533D+00

    bx = b * x
    t = x - bx

    if ( b * eps <= 0.02D+00 ) then
        c = log ( x ) + psi ( b ) + g + t
    else
        c = log ( bx ) + g + t
    end if

    tol = 5.0D+00 * eps * abs ( c )
    j = 1.0D+00
    s = 0.0D+00

    do

        j = j + 1.0D+00
        t = t * ( x - bx / j )
        aj = t / j
        s = s + aj

        if ( abs ( aj ) <= tol ) then
            exit
        end if

    end do

    apser = -a * ( c + s )

    return
    end function apser

    function bcorr ( a0, b0 )

    !*****************************************************************************80
    !
    !! BCORR evaluates DEL(A0) + DEL(B0) - DEL(A0 + B0).
    !
    !  Discussion:
    !
    !    The function DEL(A) is a remainder term that is used in the expression:
    !
    !      ln ( Gamma ( A ) ) = ( A - 0.5 ) * ln ( A )
    !        - A + 0.5 * ln ( 2 * PI ) + DEL ( A ),
    !
    !    or, in other words, DEL ( A ) is defined as:
    !
    !      DEL ( A ) = ln ( Gamma ( A ) ) - ( A - 0.5 ) * ln ( A )
    !        + A + 0.5 * ln ( 2 * PI ).
    !
    !  Author:
    !
    !    Armido DiDinato, Alfred Morris
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) A0, B0, the arguments.
    !    It is assumed that 8 <= A0 and 8 <= B0.
    !
    !    Output, real ( kind = dp ) BCORR, the value of the function.
    !
    implicit none

    real ( kind = dp ) :: a0,b0
    real ( kind = dp ) :: bcorr

    real ( kind = dp ) :: a,b,c,h,s11,s3,s5,s7,s9,t,w,x,x2

    real ( kind = dp ), parameter :: c0 =  0.833333333333333D-01
    real ( kind = dp ), parameter :: c1 = -0.277777777760991D-02
    real ( kind = dp ), parameter :: c2 =  0.793650666825390D-03
    real ( kind = dp ), parameter :: c3 = -0.595202931351870D-03
    real ( kind = dp ), parameter :: c4 =  0.837308034031215D-03
    real ( kind = dp ), parameter :: c5 = -0.165322962780713D-02

    a = min ( a0, b0 )
    b = max ( a0, b0 )

    h = a / b
    c = h / ( 1.0D+00 + h )
    x = 1.0D+00 / ( 1.0D+00 + h )
    x2 = x * x
    !
    !  Set SN = (1 - X**N)/(1 - X)
    !
    s3 = 1.0D+00 + ( x + x2 )
    s5 = 1.0D+00 + ( x + x2 * s3 )
    s7 = 1.0D+00 + ( x + x2 * s5 )
    s9 = 1.0D+00 + ( x + x2 * s7 )
    s11 = 1.0D+00 + ( x + x2 * s9 )
    !
    !  Set W = DEL(B) - DEL(A + B)
    !
    t = ( 1.0D+00 / b )**2

    w = (((( &
        c5 * s11  * t &
        + c4 * s9 ) * t &
        + c3 * s7 ) * t &
        + c2 * s5 ) * t &
        + c1 * s3 ) * t &
        + c0

    w = w * ( c / b )
    !
    !  Compute  DEL(A) + W.
    !
    t = ( 1.0D+00 / a )**2

    bcorr = ((((( &
        c5   * t &
        + c4 ) * t &
        + c3 ) * t &
        + c2 ) * t &
        + c1 ) * t &
        + c0 ) / a + w

    return
    end function bcorr

    function beta_func ( a, b )

    !*****************************************************************************80
    !
    !! BETA_FUNC evaluates the beta function.
    !
    !  Modified:
    !
    !    03 December 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) A, B, the arguments of the beta function.
    !
    !    Output, real ( kind = dp ) BETA_FUNC, the value of the beta function.
    !
    implicit none

    real ( kind = dp ) :: a,b
    real ( kind = dp ) :: beta_func

    beta_func = exp ( beta_func_log ( a, b ) )

    return
    end function beta_func

    function beta_asym ( a, b, lambda, eps )

    !*****************************************************************************80
    !
    !! BETA_ASYM computes an asymptotic expansion for IX(A,B), for large A and B.
    !
    !  Author:
    !
    !    Armido DiDinato, Alfred Morris
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) A, B, the parameters of the function.
    !    A and B should be nonnegative.  It is assumed that both A and B
    !    are greater than or equal to 15.
    !
    !    Input, real ( kind = dp ) LAMBDA, the value of ( A + B ) * Y - B.
    !    It is assumed that 0 <= LAMBDA.
    !
    !    Input, real ( kind = dp ) EPS, the tolerance.
    !
    implicit none

    real ( kind = dp ) :: a,b,lambda,eps
    real ( kind = dp ) :: beta_asym

    integer ( kind = int32 ), parameter :: num = 20

    real ( kind = dp ) a0(num+1)
    real ( kind = dp ) b0(num+1)
    real ( kind = dp ) bsum
    real ( kind = dp ) c(num+1)
    real ( kind = dp ) d(num+1)
    real ( kind = dp ) dsum

    real ( kind = dp ) f,h,h2,hn,j0,j1,r,r0,r1,s,sum1
    real ( kind = dp ) t,t0,t1,u,w,w0,z,z0,z2,zn,znm1

    integer ( kind = int32 ) i,j,m,mm1,mmj,n,np1

    real ( kind = dp ), parameter :: e0 = 1.12837916709551D+00
    real ( kind = dp ), parameter :: e1 = 0.353553390593274D+00

    beta_asym = 0.0D+00

    if ( a < b ) then
        h = a / b
        r0 = 1.0D+00 / ( 1.0D+00 + h )
        r1 = ( b - a ) / b
        w0 = 1.0D+00 / sqrt ( a * ( 1.0D+00 + h ))
    else
        h = b / a
        r0 = 1.0D+00 / ( 1.0D+00 + h )
        r1 = ( b - a ) / a
        w0 = 1.0D+00 / sqrt ( b * ( 1.0D+00 + h ))
    end if

    f = a * rlog1 ( - lambda / a ) + b * rlog1 ( lambda / b )
    t = exp ( - f )
    if ( t == 0.0D+00 ) then
        return
    end if

    z0 = sqrt ( f )
    z = 0.5D+00 * ( z0 / e1 )
    z2 = f + f

    a0(1) = ( 2.0D+00 / 3.0D+00 ) * r1
    c(1) = -0.5D+00 * a0(1)
    d(1) = -c(1)
    j0 = ( 0.5D+00 / e0 ) * error_fc ( 1, z0 )
    j1 = e1
    sum1 = j0 + d(1) * w0 * j1

    s = 1.0D+00
    h2 = h * h
    hn = 1.0D+00
    w = w0
    znm1 = z
    zn = z2

    do n = 2, num, 2

        hn = h2 * hn
        a0(n) = 2.0D+00 * r0 * ( 1.0D+00 + h * hn ) &
            / ( n + 2.0D+00 )
        np1 = n + 1
        s = s + hn
        a0(np1) = 2.0D+00 * r1 * s / ( n + 3.0D+00 )

        do i = n, np1

            r = -0.5D+00 * ( i + 1.0D+00 )
            b0(1) = r * a0(1)
            do m = 2, i
                bsum = 0.0D+00
                mm1 = m - 1
                do j = 1, mm1
                    mmj = m - j
                    bsum = bsum + ( j * r - mmj ) * a0(j) * b0(mmj)
                end do
                b0(m) = r * a0(m) + bsum / m
            end do

            c(i) = b0(i) / ( i + 1.0D+00 )

            dsum = 0.0
            do j = 1, i-1
                dsum = dsum + d(i-j) * c(j)
            end do
            d(i) = - ( dsum + c(i) )

        end do

        j0 = e1 * znm1 + ( n - 1.0D+00 ) * j0
        j1 = e1 * zn + n * j1
        znm1 = z2 * znm1
        zn = z2 * zn
        w = w0 * w
        t0 = d(n) * w * j0
        w = w0 * w
        t1 = d(np1) * w * j1
        sum1 = sum1 + ( t0 + t1 )

        if ( ( abs ( t0 ) + abs ( t1 )) <= eps * sum1 ) then
            u = exp ( - bcorr ( a, b ) )
            beta_asym = e0 * t * u * sum1
            return
        end if

    end do

    u = exp ( - bcorr ( a, b ) )
    beta_asym = e0 * t * u * sum1

    return
    end function beta_asym

    function beta_frac ( a, b, x, y, lambda, eps )

    !*****************************************************************************80
    !
    !! BETA_FRAC evaluates a continued fraction expansion for IX(A,B).
    !
    !  Author:
    !
    !    Armido DiDinato, Alfred Morris
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) A, B, the parameters of the function.
    !    A and B should be nonnegative.  It is assumed that both A and
    !    B are greater than 1.
    !
    !    Input, real ( kind = dp ) X, Y.  X is the argument of the
    !    function, and should satisy 0 <= X <= 1.  Y should equal 1 - X.
    !
    !    Input, real ( kind = dp ) LAMBDA, the value of ( A + B ) * Y - B.
    !
    !    Input, real ( kind = dp ) EPS, a tolerance.
    !
    !    Output, real ( kind = dp ) BETA_FRAC, the value of the continued
    !    fraction approximation for IX(A,B).
    !
    implicit none

    real ( kind = dp ) :: a,b,x,y,lambda,eps
    real ( kind = dp ) :: beta_frac

    real ( kind = dp ) :: alpha,an,anp1,beta,bn,bnp1,c,c0,c1
    real ( kind = dp ) :: e,n,p,r,r0,s,t,w,yp1

    beta_frac = beta_rcomp ( a, b, x, y )

    if ( beta_frac == 0.0D+00 ) then
        return
    end if

    c = 1.0D+00 + lambda
    c0 = b / a
    c1 = 1.0D+00 + 1.0D+00 / a
    yp1 = y + 1.0D+00

    n = 0.0D+00
    p = 1.0D+00
    s = a + 1.0D+00
    an = 0.0D+00
    bn = 1.0D+00
    anp1 = 1.0D+00
    bnp1 = c / c1
    r = c1 / c
    !
    !  Continued fraction calculation.
    !
    do

        n = n + 1.0D+00
        t = n / a
        w = n * ( b - n ) * x
        e = a / s
        alpha = ( p * ( p + c0 ) * e * e ) * ( w * x )
        e = ( 1.0D+00 + t ) / ( c1 + t + t )
        beta = n + w / s + e * ( c + n * yp1 )
        p = 1.0D+00 + t
        s = s +  2.0D+00
        !
        !  Update AN, BN, ANP1, and BNP1.
        !
        t = alpha * an + beta * anp1
        an = anp1
        anp1 = t
        t = alpha * bn + beta * bnp1
        bn = bnp1
        bnp1 = t

        r0 = r
        r = anp1 / bnp1

        if ( abs ( r - r0 ) <= eps * r ) then
            beta_frac = beta_frac * r
            exit
        end if
        !
        !  Rescale AN, BN, ANP1, and BNP1.
        !
        an = an / bnp1
        bn = bn / bnp1
        anp1 = r
        bnp1 = 1.0D+00

    end do

    return
    end function beta_frac

    subroutine beta_grat ( a, b, x, y, w, eps, ierr )

    !*****************************************************************************80
    !
    !! BETA_GRAT evaluates an asymptotic expansion for IX(A,B).
    !
    !  Author:
    !
    !    Armido DiDinato, Alfred Morris
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) A, B, the parameters of the function.
    !    A and B should be nonnegative.  It is assumed that 15 <= A
    !    and B <= 1, and that B is less than A.
    !
    !    Input, real ( kind = dp ) X, Y.  X is the argument of the
    !    function, and should satisy 0 <= X <= 1.  Y should equal 1 - X.
    !
    !    Input/output, real ( kind = dp ) W, a quantity to which the
    !    result of the computation is to be added on output.
    !
    !    Input, real ( kind = dp ) EPS, a tolerance.
    !
    !    Output, integer ( kind = int32 ) IERR, an error flag, which is 0 if no error
    !    was detected.
    !
    implicit none

    real ( kind = dp ) :: a,b,x,y,w,eps
    integer ( kind = int32 ) :: ierr

    real ( kind = dp ) :: bm1,bp2n,cn,coef,dj,j,l,lnx,n2,nu,p,q,r,s,sum1,t,t2,u,v,z
    real ( kind = dp ) :: c(30), d(30)

    integer ( kind = int32 ) :: i, n

    bm1 = ( b - 0.5D+00 ) - 0.5D+00
    nu = a + 0.5D+00 * bm1

    if ( y <= 0.375D+00 ) then
        lnx = alnrel ( - y )
    else
        lnx = log ( x )
    end if

    z = -nu * lnx

    if ( b * z == 0.0D+00 ) then
        ierr = 1
        return
    end if
    !
    !  Computation of the expansion.
    !
    !  Set R = EXP(-Z)*Z^B/GAMMA(B)
    !
    r = b * ( 1.0D+00 + gam1 ( b ) ) * exp ( b * log ( z ))
    r = r * exp ( a * lnx ) * exp ( 0.5D+00 * bm1 * lnx )
    u = algdiv ( b, a ) + b * log ( nu )
    u = r * exp ( - u )

    if ( u == 0.0D+00 ) then
        ierr = 1
        return
    end if

    call gamma_rat1 ( b, z, r, p, q, eps )

    v = 0.25D+00 * ( 1.0D+00 / nu )**2
    t2 = 0.25D+00 * lnx * lnx
    l = w / u
    j = q / r
    sum1 = j
    t = 1.0D+00
    cn = 1.0D+00
    n2 = 0.0D+00

    do n = 1, 30

        bp2n = b + n2
        j = ( bp2n * ( bp2n + 1.0D+00 ) * j &
            + ( z + bp2n + 1.0D+00 ) * t ) * v
        n2 = n2 +  2.0D+00
        t = t * t2
        cn = cn / ( n2 * ( n2 + 1.0D+00 ))
        c(n) = cn
        s = 0.0D+00

        coef = b - n
        do i = 1, n-1
            s = s + coef * c(i) * d(n-i)
            coef = coef + b
        end do

        d(n) = bm1 * cn + s / n
        dj = d(n) * j
        sum1 = sum1 + dj

        if ( sum1 <= 0.0D+00 ) then
            ierr = 1
            return
        end if

        if ( abs ( dj ) <= eps * ( sum1 + l ) ) then
            ierr = 0
            w = w + u * sum1
            return
        end if

    end do

    ierr = 0
    w = w + u * sum1

    return
    end subroutine beta_grat

    subroutine beta_inc ( a, b, x, y, w, w1, ierr )

    !*****************************************************************************80
    !
    !! BETA_INC evaluates the incomplete beta function IX(A,B).
    !
    !  Author:
    !
    !    Alfred Morris,
    !    Naval Surface Weapons Center,
    !    Dahlgren, Virginia.
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) A, B, the parameters of the function.
    !    A and B should be nonnegative.
    !
    !    Input, real ( kind = dp ) X, Y.  X is the argument of the
    !    function, and should satisy 0 <= X <= 1.  Y should equal 1 - X.
    !
    !    Output, real ( kind = dp ) W, W1, the values of IX(A,B) and
    !    1-IX(A,B).
    !
    !    Output, integer ( kind = int32 ) IERR, the error flag.
    !    0, no error was detected.
    !    1, A or B is negative;
    !    2, A = B = 0;
    !    3, X < 0 or 1 < X;
    !    4, Y < 0 or 1 < Y;
    !    5, X + Y /= 1;
    !    6, X = A = 0;
    !    7, Y = B = 0.
    !
    implicit none

    real ( kind = dp ) :: a,b,x,y,w,w1
    integer ( kind = int32 ) :: ierr

    real ( kind = dp ) :: a0,b0,eps,lambda,t,x0,y0,z
    integer ( kind = int32 ) :: ierr1,ind,n

    eps = epsilon ( eps )
    w = 0.0D+00
    w1 = 0.0D+00

    if ( a < 0.0D+00 .or. b < 0.0D+00 ) then
        ierr = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BETA_INC - Fatal error!'
        write ( *, '(a,i8)' ) '  IERR = ', ierr
        return
    end if

    if ( a == 0.0D+00 .and. b == 0.0D+00 ) then
        ierr = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BETA_INC - Fatal error!'
        write ( *, '(a,i8)' ) '  IERR = ', ierr
        return
    end if

    if ( x < 0.0D+00 .or. 1.0D+00 < x ) then
        ierr = 3
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BETA_INC - Fatal error!'
        write ( *, '(a,i8)' ) '  IERR = ', ierr
        return
    end if

    if ( y < 0.0D+00 .or. 1.0D+00 < y ) then
        ierr = 4
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BETA_INC - Fatal error!'
        write ( *, '(a,i8)' ) '  IERR = ', ierr
        return
    end if

    z = ( ( x + y ) - 0.5D+00 ) - 0.5D+00

    if ( 3.0D+00 * eps < abs ( z ) ) then
        ierr = 5
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BETA_INC - Fatal error!'
        write ( *, '(a,i8)' ) '  IERR = ', ierr
        return
    end if

    ierr = 0

    if ( x == 0.0D+00 ) then
        w = 0.0D+00
        w1 = 1.0D+00
        if ( a == 0.0D+00 ) then
            ierr = 6
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'BETA_INC - Fatal error!'
            write ( *, '(a,i8)' ) '  IERR = ', ierr
        end if
        return
    end if

    if ( y == 0.0D+00 ) then
        if ( b == 0.0D+00 ) then
            ierr = 7
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'BETA_INC - Fatal error!'
            write ( *, '(a,i8)' ) '  IERR = ', ierr
            return
        end if
        w = 1.0D+00
        w1 = 0.0D+00
        return
    end if

    if ( a == 0.0D+00 ) then
        w = 1.0D+00
        w1 = 0.0D+00
        return
    end if

    if ( b == 0.0D+00 ) then
        w = 0.0D+00
        w1 = 1.0D+00
        return
    end if

    eps = max ( eps, 1.0D-15 )

    if ( max ( a, b ) < 0.001D+00 * eps ) then
        go to 260
    end if

    ind = 0
    a0 = a
    b0 = b
    x0 = x
    y0 = y

    if ( 1.0D+00 < min ( a0, b0 ) ) then
        go to 40
    end if
    !
    !  Procedure for A0 <= 1 or B0 <= 1
    !
    if ( 0.5D+00 < x ) then
        ind = 1
        a0 = b
        b0 = a
        x0 = y
        y0 = x
    end if

    if ( b0 < min ( eps, eps * a0 ) ) then
        go to 90
    end if

    if ( a0 < min ( eps, eps * b0 ) .and. b0 * x0 <= 1.0D+00 ) then
        go to 100
    end if

    if ( 1.0D+00 < max ( a0, b0 ) ) then
        go to 20
    end if

    if ( min ( 0.2D+00, b0 ) <= a0 ) then
        go to 110
    end if

    if ( x0**a0 <= 0.9D+00 ) then
        go to 110
    end if

    if ( 0.3D+00 <= x0 ) then
        go to 120
    end if

    n = 20
    go to 140

20  continue

    if ( b0 <= 1.0D+00 ) then
        go to 110
    end if

    if ( 0.3D+00 <= x0 ) then
        go to 120
    end if

    if ( 0.1D+00 <= x0 ) then
        go to 30
    end if

    if ( ( x0 * b0 )**a0 <= 0.7D+00 ) then
        go to 110
    end if

30  continue

    if ( 15.0D+00 < b0 ) then
        go to 150
    end if

    n = 20
    go to 140
    !
    !  PROCEDURE for 1 < A0 and 1 < B0.
    !
40  continue

    if ( a <= b ) then
        lambda = a - ( a + b ) * x
    else
        lambda = ( a + b ) * y - b
    end if

    if ( lambda < 0.0D+00 ) then
        ind = 1
        a0 = b
        b0 = a
        x0 = y
        y0 = x
        lambda = abs ( lambda )
    end if

70  continue

    if ( b0 < 40.0D+00 .and. b0 * x0 <= 0.7D+00 ) then
        go to 110
    end if

    if ( b0 < 40.0D+00 ) then
        go to 160
    end if

    if ( b0 < a0 ) then
        go to 80
    end if

    if ( a0 <= 100.0D+00 ) then
        go to 130
    end if

    if ( 0.03D+00 * a0 < lambda ) then
        go to 130
    end if

    go to 200

80  continue

    if ( b0 <= 100.0D+00 ) then
        go to 130
    end if

    if ( 0.03D+00 * b0 < lambda ) then
        go to 130
    end if

    go to 200
    !
    !  Evaluation of the appropriate algorithm.
    !
90  continue

    w = fpser ( a0, b0, x0, eps )
    w1 = 0.5D+00 + ( 0.5D+00 - w )
    go to 250

100 continue

    w1 = apser ( a0, b0, x0, eps )
    w = 0.5D+00 + ( 0.5D+00 - w1 )
    go to 250

110 continue

    w = beta_pser ( a0, b0, x0, eps )
    w1 = 0.5D+00 + ( 0.5D+00 - w )
    go to 250

120 continue

    w1 = beta_pser ( b0, a0, y0, eps )
    w = 0.5D+00 + ( 0.5D+00 - w1 )
    go to 250

130 continue

    w = beta_frac ( a0, b0, x0, y0, lambda, 15.0D+00 * eps )
    w1 = 0.5D+00 + ( 0.5D+00 - w )
    go to 250

140 continue

    w1 = beta_up ( b0, a0, y0, x0, n, eps )
    b0 = b0 + n

150 continue

    call beta_grat ( b0, a0, y0, x0, w1, 15.0D+00 * eps, ierr1 )
    w = 0.5D+00 + ( 0.5D+00 - w1 )
    go to 250

160 continue

    n = b0
    b0 = b0 - n

    if ( b0 == 0.0D+00 ) then
        n = n - 1
        b0 = 1.0D+00
    end if

170 continue

    w = beta_up ( b0, a0, y0, x0, n, eps )

    if ( x0 <= 0.7D+00 ) then
        w = w + beta_pser ( a0, b0, x0, eps )
        w1 = 0.5D+00 + ( 0.5D+00 - w )
        go to 250
    end if

    if ( a0 <= 15.0D+00 ) then
        n = 20
        w = w + beta_up ( a0, b0, x0, y0, n, eps )
        a0 = a0 + n
    end if

190 continue

    call beta_grat ( a0, b0, x0, y0, w, 15.0D+00 * eps, ierr1 )
    w1 = 0.5D+00 + ( 0.5D+00 - w )
    go to 250

200 continue

    w = beta_asym ( a0, b0, lambda, 100.0D+00 * eps )
    w1 = 0.5D+00 + ( 0.5D+00 - w )
    go to 250
    !
    !  Termination of the procedure.
    !
250 continue

    if ( ind /= 0 ) then
        t = w
        w = w1
        w1 = t
    end if

    return
    !
    !  Procedure for A and B < 0.001 * EPS
    !
260 continue

    w = b / ( a + b )
    w1 = a / ( a + b )

    return
    end subroutine beta_inc

    subroutine beta_inc_values ( n_data, a, b, x, fx )

    !*****************************************************************************80
    !
    !! BETA_INC_VALUES returns some values of the incomplete Beta function.
    !
    !  Discussion:
    !
    !    The incomplete Beta function may be written
    !
    !      BETA_INC(A,B,X) = Integral (0 to X) T^(A-1) * (1-T)^(B-1) dT
    !                      / Integral (0 to 1) T^(A-1) * (1-T)^(B-1) dT
    !
    !    Thus,
    !
    !      BETA_INC(A,B,0.0) = 0.0
    !      BETA_INC(A,B,1.0) = 1.0
    !
    !    Note that in Mathematica, the expressions:
    !
    !      BETA[A,B]   = Integral (0 to 1) T^(A-1) * (1-T)^(B-1) dT
    !      BETA[X,A,B] = Integral (0 to X) T^(A-1) * (1-T)^(B-1) dT
    !
    !    and thus, to evaluate the incomplete Beta function requires:
    !
    !      BETA_INC(A,B,X) = BETA[X,A,B] / BETA[A,B]
    !
    !  Modified:
    !
    !    17 February 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Milton Abramowitz, Irene Stegun,
    !    Handbook of Mathematical Functions,
    !    US Department of Commerce, 1964.
    !
    !    Karl Pearson,
    !    Tables of the Incomplete Beta Function,
    !    Cambridge University Press, 1968.
    !
    !  Parameters:
    !
    !    Input/output, integer ( kind = int32 ) N_DATA.  The user sets N_DATA to 0
    !    before the first call.  On each call, the routine increments N_DATA by 1,
    !    and returns the corresponding data; when there is no more data, the
    !    output value of N_DATA will be 0 again.
    !
    !    Output, real ( kind = dp ) A, B, X, the arguments of the function.
    !
    !    Output, real ( kind = dp ) FX, the value of the function.
    !
    implicit none

    integer ( kind = int32 ) :: n_data
    real ( kind = dp ) :: a, b, x, fx

    integer ( kind = int32 ), parameter :: n_max = 30

    real ( kind = dp ), save, dimension ( n_max ) :: a_vec = (/ &
        0.5D+00,  0.5D+00,  0.5D+00,  1.0D+00, &
        1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00, &
        2.0D+00,  2.0D+00,  2.0D+00,  2.0D+00, &
        2.0D+00,  2.0D+00,  2.0D+00,  2.0D+00, &
        2.0D+00,  5.5D+00, 10.0D+00, 10.0D+00, &
        10.0D+00, 10.0D+00, 20.0D+00, 20.0D+00, &
        20.0D+00, 20.0D+00, 20.0D+00, 30.0D+00, &
        30.0D+00, 40.0D+00 /)
    real ( kind = dp ), save, dimension ( n_max ) :: b_vec = (/ &
        0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00, &
        0.5D+00,  0.5D+00,  0.5D+00,  1.0D+00, &
        2.0D+00,  2.0D+00,  2.0D+00,  2.0D+00, &
        2.0D+00,  2.0D+00,  2.0D+00,  2.0D+00, &
        2.0D+00,  5.0D+00,  0.5D+00,  5.0D+00, &
        5.0D+00, 10.0D+00,  5.0D+00, 10.0D+00, &
        10.0D+00, 20.0D+00, 20.0D+00, 10.0D+00, &
        10.0D+00, 20.0D+00 /)
    real ( kind = dp ), save, dimension ( n_max ) :: fx_vec = (/ &
        0.0637686D+00, 0.2048328D+00, 1.0000000D+00, 0.0D+00,       &
        0.0050126D+00, 0.0513167D+00, 0.2928932D+00, 0.5000000D+00, &
        0.028D+00,     0.104D+00,     0.216D+00,     0.352D+00,     &
        0.500D+00,     0.648D+00,     0.784D+00,     0.896D+00,     &
        0.972D+00,     0.4361909D+00, 0.1516409D+00, 0.0897827D+00, &
        1.0000000D+00, 0.5000000D+00, 0.4598773D+00, 0.2146816D+00, &
        0.9507365D+00, 0.5000000D+00, 0.8979414D+00, 0.2241297D+00, &
        0.7586405D+00, 0.7001783D+00 /)
    real ( kind = dp ), save, dimension ( n_max ) :: x_vec = (/ &
        0.01D+00, 0.10D+00, 1.00D+00, 0.0D+00,  &
        0.01D+00, 0.10D+00, 0.50D+00, 0.50D+00, &
        0.1D+00,  0.2D+00,  0.3D+00,  0.4D+00,  &
        0.5D+00,  0.6D+00,  0.7D+00,  0.8D+00,  &
        0.9D+00,  0.50D+00, 0.90D+00, 0.50D+00, &
        1.00D+00, 0.50D+00, 0.80D+00, 0.60D+00, &
        0.80D+00, 0.50D+00, 0.60D+00, 0.70D+00, &
        0.80D+00, 0.70D+00 /)

    if ( n_data < 0 ) then
        n_data = 0
    end if

    n_data = n_data + 1

    if ( n_max < n_data ) then
        n_data = 0
        a = 0.0D+00
        b = 0.0D+00
        x = 0.0D+00
        fx = 0.0D+00
    else
        a = a_vec(n_data)
        b = b_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
    end if

    return
    end subroutine beta_inc_values

    function beta_func_log ( a0, b0 )

    !*****************************************************************************80
    !
    !! BETA_FUNC_LOG evaluates the logarithm of the beta function.
    !
    !  Author:
    !
    !    Armido DiDinato, Alfred Morris
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) A0, B0, the parameters of the function.
    !    A0 and B0 should be nonnegative.
    !
    !    Output, real ( kind = dp ) BETA_FUNC_LOG, the value of the logarithm
    !    of the Beta function.
    !
    implicit none

    real ( kind = dp ) :: a0, b0
    real ( kind = dp ) ::  beta_func_log

    real ( kind = dp ) :: a,b,c,h,u,v,w,z
    integer ( kind = int32 ) :: i,n

    real ( kind = dp ), parameter :: e = 0.918938533204673D+00

    a = min ( a0, b0 )
    b = max ( a0, b0 )
    !
    !  8 < A.
    !
    if ( 8.0D+00 <= a ) then

        w = bcorr ( a, b )
        h = a / b
        c = h / ( 1.0D+00 + h )
        u = - ( a - 0.5D+00 ) * log ( c )
        v = b * alnrel ( h )

        if ( v < u ) then
            beta_func_log = ((( -0.5D+00 * log ( b ) + e ) + w ) - v ) - u
        else
            beta_func_log = ((( -0.5D+00 * log ( b ) + e ) + w ) - u ) - v
        end if

        return
    end if
    !
    !  Procedure when A < 1
    !
    if ( a < 1.0D+00 ) then

        if ( b < 8.0D+00 ) then
            beta_func_log = gamma_log ( a ) + ( gamma_log ( b ) - gamma_log ( a + b ) )
        else
            beta_func_log = gamma_log ( a ) + algdiv ( a, b )
        end if

        return

    end if
    !
    !  Procedure when 1 <= A < 8
    !
    if ( 2.0D+00 < a ) then
        go to 40
    end if

    if ( b <= 2.0D+00 ) then
        beta_func_log = gamma_log ( a ) + gamma_log ( b ) - gsumln ( a, b )
        return
    end if

    w = 0.0D+00

    if ( b < 8.0D+00 ) then
        go to 60
    end if

    beta_func_log = gamma_log ( a ) + algdiv ( a, b )
    return

40  continue
    !
    !  Reduction of A when 1000 < B.
    !
    if ( 1000.0D+00 < b ) then

        n = a - 1.0D+00
        w = 1.0D+00
        do i = 1, n
            a = a - 1.0D+00
            w = w * ( a / ( 1.0D+00 + a / b ))
        end do

        beta_func_log = ( log ( w ) - n * log ( b ) ) &
            + ( gamma_log ( a ) + algdiv ( a, b ) )

        return
    end if

    n = a - 1.0D+00
    w = 1.0D+00
    do i = 1, n
        a = a - 1.0D+00
        h = a / b
        w = w * ( h / ( 1.0D+00 + h ) )
    end do
    w = log ( w )

    if ( 8.0D+00 <= b ) then
        beta_func_log = w + gamma_log ( a ) + algdiv ( a, b )
        return
    end if
    !
    !  Reduction of B when B < 8.
    !
60  continue

    n = b - 1.0D+00
    z = 1.0D+00
    do i = 1, n
        b = b - 1.0D+00
        z = z * ( b / ( a + b ))
    end do

    beta_func_log = w + log ( z ) + ( gamma_log ( a ) + ( gamma_log ( b ) &
        - gsumln ( a, b ) ) )

    return
    end function beta_func_log

    function beta_pser ( a, b, x, eps )

    !*****************************************************************************80
    !
    !! BETA_PSER uses a power series expansion to evaluate IX(A,B)(X).
    !
    !  Discussion:
    !
    !    BETA_PSER is used when B <= 1 or B*X <= 0.7.
    !
    !  Author:
    !
    !    Armido DiDinato, Alfred Morris
    !
    !    Edited in October 2018 by Patrick Macnamara to prevent routine from getting
    !    stuck in an infinite loop if one of the inputs is NaN.
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) A, B, the parameters.
    !
    !    Input, real ( kind = dp ) X, the point where the function
    !    is to be evaluated.
    !
    !    Input, real ( kind = dp ) EPS, the tolerance.
    !
    !    Output, real ( kind = dp ) BETA_PSER, the approximate value of IX(A,B)(X).
    !
    implicit none

    real ( kind = dp ) :: a,b,x,eps
    real ( kind = dp ) :: beta_pser

    real ( kind = dp ) :: a0,b0,apb,c,n,sum1,t,tol,u,w,z

    integer ( kind = int32 ) :: i,m

    beta_pser = 0.0D+00

    if ( x == 0.0D+00 ) then
        return
    end if
    !
    !  Compute the factor X**A/(A*BETA(A,B))
    !
    a0 = min ( a, b )

    if ( 1.0D+00 <= a0 ) then

        z = a * log ( x ) - beta_func_log ( a, b )
        beta_pser = exp ( z ) / a

    else

        b0 = max ( a, b )

        if ( b0 <= 1.0D+00 ) then

            beta_pser = x**a
            if ( beta_pser == 0.0D+00 ) then
                return
            end if

            apb = a + b

            if ( apb <= 1.0D+00 ) then
                z = 1.0D+00 + gam1 ( apb )
            else
                u = a + b - 1.0D+00
                z = ( 1.0D+00 + gam1 ( u ) ) / apb
            end if

            c = ( 1.0D+00 + gam1 ( a ) ) &
                * ( 1.0D+00 + gam1 ( b ) ) / z
            beta_pser = beta_pser * c * ( b / apb )

        else if ( b0 < 8.0D+00 ) then

            u = gamma_ln1 ( a0 )
            m = b0 - 1.0D+00

            c = 1.0D+00
            do i = 1, m
                b0 = b0 - 1.0D+00
                c = c * ( b0 / ( a0 + b0 ))
            end do

            u = log ( c ) + u
            z = a * log ( x ) - u
            b0 = b0 - 1.0D+00
            apb = a0 + b0

            if ( apb <= 1.0D+00 ) then
                t = 1.0D+00 + gam1 ( apb )
            else
                u = a0 + b0 - 1.0D+00
                t = ( 1.0D+00 + gam1 ( u ) ) / apb
            end if

            beta_pser = exp ( z ) * ( a0 / a ) &
                * ( 1.0D+00 + gam1 ( b0 )) / t

        else if ( 8.0D+00 <= b0 ) then

            u = gamma_ln1 ( a0 ) + algdiv ( a0, b0 )
            z = a * log ( x ) - u
            beta_pser = ( a0 / a ) * exp ( z )

        end if

    end if

    if ( beta_pser == 0.0D+00 .or. a <= 0.1D+00 * eps ) then
        return
    end if
    !
    !  Compute the series.
    !
    sum1 = 0.0D+00
    n = 0.0D+00
    c = 1.0D+00
    tol = eps / a

    do

        n = n + 1.0D+00
        c = c * ( 0.5D+00 + ( 0.5D+00 - b / n ) ) * x
        w = c / ( a + n )
        sum1 = sum1 + w

        if ( abs ( w ) <= tol ) then
            exit
        end if

        ! Added by Patrick Macnamara to prevent routine getting stuck in an
        ! infinite loop if w is NaN
        ! exit if w is NaN
        if ( w /= w ) then
            exit
        end if
    end do

    beta_pser = beta_pser * ( 1.0D+00 + a * sum1 )

    return
    end function beta_pser

    function beta_rcomp ( a, b, x, y )

    !*****************************************************************************80
    !
    !! BETA_RCOMP evaluates X^A * Y^B / Beta(A,B).
    !
    !  Author:
    !
    !    Armido DiDinato, Alfred Morris
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) A, B, the parameters of the Beta function.
    !    A and B should be nonnegative.
    !
    !    Input, real ( kind = dp ) X, Y, define the numerator of the fraction.
    !
    !    Output, real ( kind = dp ) BETA_RCOMP, the value of X^A * Y^B / Beta(A,B).
    !
    implicit none

    real ( kind = dp ) :: a,b,x,y
    real ( kind = dp ) :: beta_rcomp

    real ( kind = dp ) :: a0,b0,apb,c,e,h,lambda,lnx,lny,t,u,v,x0,y0,z
    integer ( kind = int32 ) :: i,n

    real ( kind = dp ), parameter :: const = 0.398942280401433D+00

    beta_rcomp = 0.0D+00
    if ( x == 0.0D+00 .or. y == 0.0D+00 ) then
        return
    end if

    a0 = min ( a, b )

    if ( a0 < 8.0D+00 ) then

        if ( x <= 0.375D+00 ) then
            lnx = log ( x )
            lny = alnrel ( - x )
        else if ( y <= 0.375D+00 ) then
            lnx = alnrel ( - y )
            lny = log ( y )
        else
            lnx = log ( x )
            lny = log ( y )
        end if

        z = a * lnx + b * lny

        if ( 1.0D+00 <= a0 ) then
            z = z - beta_func_log ( a, b )
            beta_rcomp = exp ( z )
            return
        end if
        !
        !  Procedure for A < 1 or B < 1
        !
        b0 = max ( a, b )

        if ( b0 <= 1.0D+00 ) then

            beta_rcomp = exp ( z )
            if ( beta_rcomp == 0.0D+00 ) then
                return
            end if

            apb = a + b

            if ( apb <= 1.0D+00 ) then
                z = 1.0D+00 + gam1 ( apb )
            else
                u = a + b - 1.0D+00
                z = ( 1.0D+00 + gam1 ( u ) ) / apb
            end if

            c = ( 1.0D+00 + gam1 ( a ) ) &
                * ( 1.0D+00 + gam1 ( b ) ) / z
            beta_rcomp = beta_rcomp * ( a0 * c ) &
                / ( 1.0D+00 + a0 / b0 )

        else if ( b0 < 8.0D+00 ) then

            u = gamma_ln1 ( a0 )
            n = b0 - 1.0D+00

            c = 1.0D+00
            do i = 1, n
                b0 = b0 - 1.0D+00
                c = c * ( b0 / ( a0 + b0 ))
            end do
            u = log ( c ) + u

            z = z - u
            b0 = b0 - 1.0D+00
            apb = a0 + b0

            if ( apb <= 1.0D+00 ) then
                t = 1.0D+00 + gam1 ( apb )
            else
                u = a0 + b0 - 1.0D+00
                t = ( 1.0D+00 + gam1 ( u ) ) / apb
            end if

            beta_rcomp = a0 * exp ( z ) * ( 1.0D+00 + gam1 ( b0 ) ) / t

        else if ( 8.0D+00 <= b0 ) then

            u = gamma_ln1 ( a0 ) + algdiv ( a0, b0 )
            beta_rcomp = a0 * exp ( z - u )

        end if

    else

        if ( a <= b ) then
            h = a / b
            x0 = h / ( 1.0D+00 + h )
            y0 = 1.0D+00 / (  1.0D+00 + h )
            lambda = a - ( a + b ) * x
        else
            h = b / a
            x0 = 1.0D+00 / ( 1.0D+00 + h )
            y0 = h / ( 1.0D+00 + h )
            lambda = ( a + b ) * y - b
        end if

        e = -lambda / a

        if ( abs ( e ) <= 0.6D+00 ) then
            u = rlog1 ( e )
        else
            u = e - log ( x / x0 )
        end if

        e = lambda / b

        if ( abs ( e ) <= 0.6D+00 ) then
            v = rlog1 ( e )
        else
            v = e - log ( y / y0 )
        end if

        z = exp ( - ( a * u + b * v ) )
        beta_rcomp = const * sqrt ( b * x0 ) * z * exp ( - bcorr ( a, b ))

    end if

    return
    end function beta_rcomp

    function beta_rcomp1 ( mu, a, b, x, y )

    !*****************************************************************************80
    !
    !! BETA_RCOMP1 evaluates exp(MU) * X^A * Y^B / Beta(A,B).
    !
    !  Author:
    !
    !    Armido DiDinato, Alfred Morris
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = int32 ) MU, ?
    !
    !    Input, real ( kind = dp ) A, B, the parameters of the Beta function.
    !    A and B should be nonnegative.
    !
    !    Input, real ( kind = dp ) X, Y, quantities whose powers form part of
    !    the expression.
    !
    !    Output, real ( kind = dp ) BETA_RCOMP1, the value of
    !    exp(MU) * X**A * Y**B / Beta(A,B).
    !
    implicit none

    integer ( kind = int32 ) :: mu
    real ( kind = dp ) :: a,b,x,y
    real ( kind = dp ) :: beta_rcomp1

    real ( kind = dp ) :: a0,b0,apb,c,e,h,lambda,lnx,lny
    real ( kind = dp ) :: t,u,v,x0,y0,z
    integer ( kind = int32 ) :: i,n

    real ( kind = dp ), parameter :: const = 0.398942280401433D+00

    a0 = min ( a, b )
    !
    !  Procedure for 8 <= A and 8 <= B.
    !
    if ( 8.0D+00 <= a0 ) then

        if ( a <= b ) then
            h = a / b
            x0 = h / ( 1.0D+00 + h )
            y0 = 1.0D+00 / ( 1.0D+00 + h )
            lambda = a - ( a + b ) * x
        else
            h = b / a
            x0 = 1.0D+00 / ( 1.0D+00 + h )
            y0 = h / ( 1.0D+00 + h )
            lambda = ( a + b ) * y - b
        end if

        e = -lambda / a

        if ( abs ( e ) <= 0.6D+00 ) then
            u = rlog1 ( e )
        else
            u = e - log ( x / x0 )
        end if

        e = lambda / b

        if ( abs ( e ) <= 0.6D+00 ) then
            v = rlog1 ( e )
        else
            v = e - log ( y / y0 )
        end if

        z = esum ( mu, - ( a * u + b * v ))
        beta_rcomp1 = const * sqrt ( b * x0 ) * z * exp ( - bcorr ( a, b ) )
        !
        !  Procedure for A < 8 or B < 8.
        !
    else

        if ( x <= 0.375D+00 ) then
            lnx = log ( x )
            lny = alnrel ( - x )
        else if ( y <= 0.375D+00 ) then
            lnx = alnrel ( - y )
            lny = log ( y )
        else
            lnx = log ( x )
            lny = log ( y )
        end if

        z = a * lnx + b * lny

        if ( 1.0D+00 <= a0 ) then
            z = z - beta_func_log ( a, b )
            beta_rcomp1 = esum ( mu, z )
            return
        end if
        !
        !  Procedure for A < 1 or B < 1.
        !
        b0 = max ( a, b )

        if ( 8.0D+00 <= b0 ) then
            u = gamma_ln1 ( a0 ) + algdiv ( a0, b0 )
            beta_rcomp1 = a0 * esum ( mu, z-u )
            return
        end if

        if ( 1.0D+00 < b0 ) then
            !
            !  Algorithm for 1 < B0 < 8
            !
            u = gamma_ln1 ( a0 )
            n = b0 - 1.0D+00

            c = 1.0D+00
            do i = 1, n
                b0 = b0 - 1.0D+00
                c = c * ( b0 / ( a0 + b0 ) )
            end do
            u = log ( c ) + u

            z = z - u
            b0 = b0 - 1.0D+00
            apb = a0 + b0

            if ( apb <= 1.0D+00 ) then
                t = 1.0D+00 + gam1 ( apb )
            else
                u = a0 + b0 - 1.0D+00
                t = ( 1.0D+00 + gam1 ( u ) ) / apb
            end if

            beta_rcomp1 = a0 * esum ( mu, z ) &
                * ( 1.0D+00 + gam1 ( b0 ) ) / t
            !
            !  Algorithm for B0 <= 1
            !
        else

            beta_rcomp1 = esum ( mu, z )
            if ( beta_rcomp1 == 0.0D+00 ) then
                return
            end if

            apb = a + b

            if ( apb <= 1.0D+00 ) then
                z = 1.0D+00 + gam1 ( apb )
            else
                u = real ( a, kind = dp ) + real ( b, kind = dp ) - 1.0D+00
                z = ( 1.0D+00 + gam1 ( u )) / apb
            end if

            c = ( 1.0D+00 + gam1 ( a ) ) &
                * ( 1.0D+00 + gam1 ( b ) ) / z
            beta_rcomp1 = beta_rcomp1 * ( a0 * c ) / ( 1.0D+00 + a0 / b0 )

        end if

    end if

    return
    end function beta_rcomp1

    function beta_up ( a, b, x, y, n, eps )

    !*****************************************************************************80
    !
    !! BETA_UP evaluates IX(A,B) - IX(A+N,B) where N is a positive integer.
    !
    !  Author:
    !
    !    Armido DiDinato, Alfred Morris
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) A, B, the parameters of the function.
    !    A and B should be nonnegative.
    !
    !    Input, real ( kind = dp ) X, Y, ?
    !
    !    Input, integer ( kind = int32 ) N, the increment to the first argument of IX.
    !
    !    Input, real ( kind = dp ) EPS, the tolerance.
    !
    !    Output, real ( kind = dp ) BETA_UP, the value of IX(A,B) - IX(A+N,B).
    !
    implicit none

    real ( kind = dp ) :: a,b,x,y
    integer ( kind = int32 ) :: n
    real ( kind = dp ) :: eps
    real ( kind = dp ) :: beta_up

    real ( kind = dp ) :: ap1,apb,d,l,r,t,w
    integer ( kind = int32 )  :: i,k,mu

    !
    !  Obtain the scaling factor EXP(-MU) AND
    !  EXP(MU)*(X**A*Y**B/BETA(A,B))/A
    !
    apb = a + b
    ap1 = a + 1.0D+00
    mu = 0
    d = 1.0D+00

    if ( n /= 1 ) then

        if ( 1.0D+00 <= a ) then

            if ( 1.1D+00 * ap1 <= apb ) then
                mu = abs ( exparg ( 1 ) )
                k = exparg ( 0 )
                if ( k < mu ) then
                    mu = k
                end if
                t = mu
                d = exp ( - t )
            end if

        end if

    end if

    beta_up = beta_rcomp1 ( mu, a, b, x, y ) / a

    if ( n == 1 .or. beta_up == 0.0D+00 ) then
        return
    end if

    w = d
    !
    !  Let K be the index of the maximum term.
    !
    k = 0

    if ( 1.0D+00 < b ) then

        if ( y <= 0.0001D+00 ) then

            k = n - 1

        else

            r = ( b - 1.0D+00 ) * x / y - a

            if ( 1.0D+00 <= r ) then
                k = n - 1
                t = n - 1
                if ( r < t ) then
                    k = r
                end if
            end if

        end if
        !
        !  Add the increasing terms of the series.
        !
        do i = 1, k
            l = i - 1
            d = ( ( apb + l ) / ( ap1 + l ) ) * x * d
            w = w + d
        end do

    end if
    !
    !  Add the remaining terms of the series.
    !
    do i = k+1, n-1
        l = i - 1
        d = ( ( apb + l ) / ( ap1 + l ) ) * x * d
        w = w + d
        if ( d <= eps * w ) then
            beta_up = beta_up * w
            return
        end if
    end do

    beta_up = beta_up * w

    return
    end function beta_up

    subroutine binomial_cdf_values ( n_data, a, b, x, fx )

    !*****************************************************************************80
    !
    !! BINOMIAL_CDF_VALUES returns some values of the binomial CDF.
    !
    !  Discussion:
    !
    !    CDF(X)(A,B) is the probability of at most X successes in A trials,
    !    given that the probability of success on a single trial is B.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    27 May 2001
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Milton Abramowitz, Irene Stegun,
    !    Handbook of Mathematical Functions,
    !    US Department of Commerce, 1964.
    !
    !    Daniel Zwillinger,
    !    CRC Standard Mathematical Tables and Formulae,
    !    30th Edition, CRC Press, 1996, pages 651-652.
    !
    !  Parameters:
    !
    !    Input/output, integer ( kind = int32 ) N_DATA.  The user sets N_DATA to 0
    !    before the first call.  On each call, the routine increments N_DATA by 1,
    !    and returns the corresponding data; when there is no more data, the
    !    output value of N_DATA will be 0 again.
    !
    !    Output, integer ( kind = int32 ) A, real ( kind = dp ) B, integer X, the
    !    arguments of the function.
    !
    !    Output, real ( kind = dp ) FX, the value of the function.
    !
    implicit none

    integer ( kind = int32 ) :: n_data
    integer ( kind = int32 ) :: a
    real ( kind = dp ) :: b
    integer ( kind = int32 ) :: x
    real ( kind = dp ) :: fx

    integer ( kind = int32 ), parameter :: n_max = 17
    integer ( kind = int32 ), save, dimension ( n_max ) :: a_vec = (/ &
        2,  2,  2,  2, &
        2,  4,  4,  4, &
        4, 10, 10, 10, &
        10, 10, 10, 10, &
        10 /)
    real ( kind = dp ), save, dimension ( n_max ) :: b_vec = (/ &
        0.05D+00, 0.05D+00, 0.05D+00, 0.50D+00, &
        0.50D+00, 0.25D+00, 0.25D+00, 0.25D+00, &
        0.25D+00, 0.05D+00, 0.10D+00, 0.15D+00, &
        0.20D+00, 0.25D+00, 0.30D+00, 0.40D+00, &
        0.50D+00 /)
    real ( kind = dp ), save, dimension ( n_max ) :: fx_vec = (/ &
        0.9025D+00, 0.9975D+00, 1.0000D+00, 0.2500D+00, &
        0.7500D+00, 0.3164D+00, 0.7383D+00, 0.9492D+00, &
        0.9961D+00, 0.9999D+00, 0.9984D+00, 0.9901D+00, &
        0.9672D+00, 0.9219D+00, 0.8497D+00, 0.6331D+00, &
        0.3770D+00 /)
    integer ( kind = int32 ), save, dimension ( n_max ) :: x_vec = (/ &
        0, 1, 2, 0, &
        1, 0, 1, 2, &
        3, 4, 4, 4, &
        4, 4, 4, 4, &
        4 /)

    if ( n_data < 0 ) then
        n_data = 0
    end if

    n_data = n_data + 1

    if ( n_max < n_data ) then
        n_data = 0
        a = 0
        b = 0.0D+00
        x = 0
        fx = 0.0D+00
    else
        a = a_vec(n_data)
        b = b_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
    end if

    return
    end subroutine binomial_cdf_values

    subroutine chi_noncentral_cdf_values ( n_data, x, lambda, df, cdf )

    !*****************************************************************************80
    !
    !! CHI_NONCENTRAL_CDF_VALUES returns values of the noncentral chi CDF.
    !
    !  Discussion:
    !
    !    The CDF of the noncentral chi square distribution can be evaluated
    !    within Mathematica by commands such as:
    !
    !      Needs["Statistics`ContinuousDistributions`"]
    !      CDF [ NoncentralChiSquareDistribution [ DF, LAMBDA ], X ]
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    12 June 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Stephen Wolfram,
    !    The Mathematica Book,
    !    Fourth Edition,
    !    Wolfram Media / Cambridge University Press, 1999.
    !
    !  Parameters:
    !
    !    Input/output, integer ( kind = int32 ) N_DATA.  The user sets N_DATA to 0
    !    before the first call.  On each call, the routine increments N_DATA by 1,
    !    and returns the corresponding data; when there is no more data, the
    !    output value of N_DATA will be 0 again.
    !
    !    Output, real ( kind = dp ) X, the argument of the function.
    !
    !    Output, real ( kind = dp ) LAMBDA, the noncentrality parameter.
    !
    !    Output, integer ( kind = int32 ) DF, the number of degrees of freedom.
    !
    !    Output, real ( kind = dp ) CDF, the noncentral chi CDF.
    !
    implicit none

    integer ( kind = int32 ) :: n_data
    real ( kind = dp ) :: x,lambda
    integer ( kind = int32 ) :: df
    real ( kind = dp ) :: cdf

    integer ( kind = int32 ), parameter :: n_max = 27
    real ( kind = dp ), save, dimension ( n_max ) :: cdf_vec = (/ &
        0.839944D+00, 0.695906D+00, 0.535088D+00, &
        0.764784D+00, 0.620644D+00, 0.469167D+00, &
        0.307088D+00, 0.220382D+00, 0.150025D+00, &
        0.307116D-02, 0.176398D-02, 0.981679D-03, &
        0.165175D-01, 0.202342D-03, 0.498448D-06, &
        0.151325D-01, 0.209041D-02, 0.246502D-03, &
        0.263684D-01, 0.185798D-01, 0.130574D-01, &
        0.583804D-01, 0.424978D-01, 0.308214D-01, &
        0.105788D+00, 0.794084D-01, 0.593201D-01 /)
    integer ( kind = int32 ), save, dimension ( n_max ) :: df_vec = (/ &
        1,   2,   3, &
        1,   2,   3, &
        1,   2,   3, &
        1,   2,   3, &
        60,  80, 100, &
        1,   2,   3, &
        10,  10,  10, &
        10,  10,  10, &
        10,  10,  10 /)
    real ( kind = dp ), save, dimension ( n_max ) :: lambda_vec = (/ &
        0.5D+00,  0.5D+00,  0.5D+00, &
        1.0D+00,  1.0D+00,  1.0D+00, &
        5.0D+00,  5.0D+00,  5.0D+00, &
        20.0D+00, 20.0D+00, 20.0D+00, &
        30.0D+00, 30.0D+00, 30.0D+00, &
        5.0D+00,  5.0D+00,  5.0D+00, &
        2.0D+00,  3.0D+00,  4.0D+00, &
        2.0D+00,  3.0D+00,  4.0D+00, &
        2.0D+00,  3.0D+00,  4.0D+00 /)
    real ( kind = dp ), save, dimension ( n_max ) :: x_vec = (/ &
        3.000D+00,  3.000D+00,  3.000D+00, &
        3.000D+00,  3.000D+00,  3.000D+00, &
        3.000D+00,  3.000D+00,  3.000D+00, &
        3.000D+00,  3.000D+00,  3.000D+00, &
        60.000D+00, 60.000D+00, 60.000D+00, &
        0.050D+00,  0.050D+00,  0.050D+00, &
        4.000D+00,  4.000D+00,  4.000D+00, &
        5.000D+00,  5.000D+00,  5.000D+00, &
        6.000D+00,  6.000D+00,  6.000D+00 /)

    if ( n_data < 0 ) then
        n_data = 0
    end if

    n_data = n_data + 1

    if ( n_max < n_data ) then
        n_data = 0
        x = 0.0D+00
        lambda = 0.0D+00
        df = 0
        cdf = 0.0D+00
    else
        x = x_vec(n_data)
        lambda = lambda_vec(n_data)
        df = df_vec(n_data)
        cdf = cdf_vec(n_data)
    end if

    return
    end subroutine chi_noncentral_cdf_values

    subroutine chi_square_cdf_values ( n_data, a, x, fx )

    !*****************************************************************************80
    !
    !! CHI_SQUARE_CDF_VALUES returns some values of the Chi-Square CDF.
    !
    !  Discussion:
    !
    !    The value of CHI_CDF ( DF, X ) can be evaluated in Mathematica by
    !    commands like:
    !
    !      Needs["Statistics`ContinuousDistributions`"]
    !      CDF[ChiSquareDistribution[DF], X ]
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    11 June 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Milton Abramowitz, Irene Stegun,
    !    Handbook of Mathematical Functions,
    !    US Department of Commerce, 1964.
    !
    !    Stephen Wolfram,
    !    The Mathematica Book,
    !    Fourth Edition,
    !    Wolfram Media / Cambridge University Press, 1999.
    !
    !  Parameters:
    !
    !    Input/output, integer ( kind = int32 ) N_DATA.  The user sets N_DATA to 0
    !    before the first call.  On each call, the routine increments N_DATA by 1,
    !    and returns the corresponding data; when there is no more data, the
    !    output value of N_DATA will be 0 again.
    !
    !    Output, integer ( kind = int32 ) A, real ( kind = dp ) X, the arguments of
    !    the function.
    !
    !    Output, real ( kind = dp ) FX, the value of the function.
    !
    implicit none

    integer ( kind = int32 ) :: n_data, a
    real ( kind = dp ) :: x, fx

    integer ( kind = int32 ), parameter :: n_max = 21
    integer ( kind = int32 ), save, dimension ( n_max ) :: a_vec = (/ &
        1,  2,  1,  2, &
        1,  2,  3,  4, &
        1,  2,  3,  4, &
        5,  3,  3,  3, &
        3,  3, 10, 10, &
        10 /)
    real ( kind = dp ), save, dimension ( n_max ) :: fx_vec = (/ &
        0.0796557D+00, 0.00498752D+00, 0.112463D+00,    0.00995017D+00, &
        0.472911D+00,  0.181269D+00,   0.0597575D+00,   0.0175231D+00, &
        0.682689D+00,  0.393469D+00,   0.198748D+00,    0.090204D+00, &
        0.0374342D+00, 0.427593D+00,   0.608375D+00,    0.738536D+00, &
        0.828203D+00,  0.88839D+00,    0.000172116D+00, 0.00365985D+00, &
        0.0185759D+00 /)
    real ( kind = dp ), save, dimension ( n_max ) :: x_vec = (/ &
        0.01D+00, 0.01D+00, 0.02D+00, 0.02D+00, &
        0.40D+00, 0.40D+00, 0.40D+00, 0.40D+00, &
        1.00D+00, 1.00D+00, 1.00D+00, 1.00D+00, &
        1.00D+00, 2.00D+00, 3.00D+00, 4.00D+00, &
        5.00D+00, 6.00D+00, 1.00D+00, 2.00D+00, &
        3.00D+00 /)

    if ( n_data < 0 ) then
        n_data = 0
    end if

    n_data = n_data + 1

    if ( n_max < n_data ) then
        n_data = 0
        a = 0
        x = 0.0D+00
        fx = 0.0D+00
    else
        a = a_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
    end if

    return
    end subroutine chi_square_cdf_values

    subroutine cumbet ( x, y, a, b, cum, ccum )

    !*****************************************************************************80
    !
    !! CUMBET evaluates the cumulative incomplete beta distribution.
    !
    !  Discussion:
    !
    !    This routine calculates the CDF to X of the incomplete beta distribution
    !    with parameters A and B.  This is the integral from 0 to x
    !    of (1/B(a,b))*f(t)) where f(t) = t**(a-1) * (1-t)**(b-1)
    !
    !  Author:
    !
    !    Barry Brown, James Lovato, Kathy Russell
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) X, the upper limit of integration.
    !
    !    Input, real ( kind = dp ) Y, the value of 1-X.
    !
    !    Input, real ( kind = dp ) A, B, the parameters of the distribution.
    !
    !    Output, real ( kind = dp ) CUM, CCUM, the values of the cumulative
    !    density function and complementary cumulative density function.
    !
    implicit none

    real ( kind = dp ) :: x,y,a,b,cum,ccum

    integer ( kind = int32 ) :: ierr

    if ( x <= 0.0D+00 ) then

        cum = 0.0
        ccum = 1.0D+00

    else if ( y <= 0.0D+00 ) then

        cum = 1.0D+00
        ccum = 0.0

    else if ( x /= x .or. y /= y ) then
        ! added by Patrick Macnamara to protect against
        ! scenario in which x or y is NaN
        cum = 0.0d0/0.0d0
        ccum = 0.0d0/0.0d0
    else

        call beta_inc ( a, b, x, y, cum, ccum, ierr )

    end if

    return
    end subroutine cumbet

    subroutine cumbin ( s, xn, pr, ompr, cum, ccum )

    !*****************************************************************************80
    !
    !! CUMBIN evaluates the cumulative binomial distribution.
    !
    !  Discussion:
    !
    !    This routine returns the probability of 0 to S successes in XN binomial
    !    trials, each of which has a probability of success, PR.
    !
    !  Author:
    !
    !    Barry Brown, James Lovato, Kathy Russell
    !
    !  Reference:
    !
    !    Milton Abramowitz, Irene Stegun,
    !    Handbook of Mathematical Functions
    !    1966, Formula 26.5.24.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) S, the upper limit of summation.
    !
    !    Input, real ( kind = dp ) XN, the number of trials.
    !
    !    Input, real ( kind = dp ) PR, the probability of success in one trial.
    !
    !    Input, real ( kind = dp ) OMPR, equals ( 1 - PR ).
    !
    !    Output, real ( kind = dp ) CUM, the cumulative binomial distribution.
    !
    !    Output, real ( kind = dp ) CCUM, the complement of the cumulative
    !    binomial distribution.
    !
    implicit none

    real ( kind = dp ) :: s,xn,pr,ompr,cum,ccum

    if ( s < xn ) then

        call cumbet ( pr, ompr, s + 1.0D+00, xn - s, ccum, cum )

    else

        cum = 1.0D+00
        ccum = 0.0D+00

    end if

    return
    end subroutine cumbin

    subroutine cumchi ( x, df, cum, ccum )

    !*****************************************************************************80
    !
    !! CUMCHI evaluates the cumulative chi-square distribution.
    !
    !  Author:
    !
    !    Barry Brown, James Lovato, Kathy Russell
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) X, the upper limit of integration.
    !
    !    Input, real ( kind = dp ) DF, the degrees of freedom of the
    !    chi-square distribution.
    !
    !    Output, real ( kind = dp ) CUM, the cumulative chi-square distribution.
    !
    !    Output, real ( kind = dp ) CCUM, the complement of the cumulative
    !    chi-square distribution.
    !
    implicit none

    real ( kind = dp ) :: x,df,cum,ccum

    real ( kind = dp ) :: a,xx

    a = df * 0.5D+00
    xx = x * 0.5D+00

    call cumgam ( xx, a, cum, ccum )

    return
    end subroutine cumchi

    subroutine cumchn ( x, df, pnonc, cum, ccum )

    !*****************************************************************************80
    !
    !! CUMCHN evaluates the cumulative noncentral chi-square distribution.
    !
    !  Discussion:
    !
    !    This routine calculates the cumulative noncentral chi-square
    !    distribution, i.e., the probability that a random variable
    !    which follows the noncentral chi-square distribution, with
    !    noncentrality parameter PNONC and continuous degrees of
    !    freedom DF, is less than or equal to X.
    !
    !  Author:
    !
    !    Barry Brown, James Lovato, Kathy Russell
    !
    !  Reference:
    !
    !    Milton Abramowitz, Irene Stegun,
    !    Handbook of Mathematical Functions
    !    1966, Formula 26.4.25.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) X, the upper limit of integration.
    !
    !    Input, real ( kind = dp ) DF, the number of degrees of freedom.
    !
    !    Input, real ( kind = dp ) PNONC, the noncentrality parameter of
    !    the noncentral chi-square distribution.
    !
    !    Output, real ( kind = dp ) CUM, CCUM, the CDF and complementary
    !    CDF of the noncentral chi-square distribution.
    !
    !  Local Parameters:
    !
    !    Local, real ( kind = dp ) EPS, the convergence criterion.  The sum
    !    stops when a term is less than EPS * SUM.
    !
    !    Local, integer NTIRED, the maximum number of terms to be evaluated
    !    in each sum.
    !
    !    Local, logical QCONV, is TRUE if convergence was achieved, that is,
    !    the program did not stop on NTIRED criterion.
    !
    implicit none

    real ( kind = dp ) :: x,df,pnonc,cum,ccum

    real ( kind = dp ) :: adj,centaj,centwt,chid2,dfd2,lcntaj,lcntwt
    real ( kind = dp ) :: lfact,pcent,pterm,sum1,sumadj,term,wt,xnonc
    integer ( kind = int32 ) :: i,icent,iterb,iterf

    real ( kind = dp ), parameter :: eps = 0.00001D+00
    integer ( kind = int32 ), parameter :: ntired = 1000

    !qsmall ( xx ) = sum1 < 1.0D-20 .or. xx < eps * sum1
    !dg(i) = df +  2.0D+00  * real ( i, kind = dp )

    if ( x <= 0.0D+00 ) then
        cum = 0.0D+00
        ccum = 1.0D+00
        return
    end if
    !
    !  When the noncentrality parameter is (essentially) zero,
    !  use cumulative chi-square distribution
    !
    if ( pnonc <= 1.0D-10 ) then
        call cumchi ( x, df, cum, ccum )
        return
    end if

    xnonc = pnonc /  2.0D+00
    !
    !  The following code calculates the weight, chi-square, and
    !  adjustment term for the central term in the infinite series.
    !  The central term is the one in which the poisson weight is
    !  greatest.  The adjustment term is the amount that must
    !  be subtracted from the chi-square to move up two degrees
    !  of freedom.
    !
    icent = int ( xnonc )
    if ( icent == 0 ) then
        icent = 1
    end if

    chid2 = x /  2.0D+00
    !
    !  Calculate central weight term.
    !
    lfact = gamma_log ( real ( icent + 1, kind = dp ) )
    lcntwt = - xnonc + icent * log ( xnonc ) - lfact
    centwt = exp ( lcntwt )
    !
    !  Calculate central chi-square.
    !
    !call cumchi ( x, dg(icent), pcent, ccum )
    call cumchi ( x, df +  2.0D+00  * real ( icent, kind = dp ), pcent, ccum )
    !
    !  Calculate central adjustment term.
    !
    !dfd2 = dg(icent) /  2.0D+00
    dfd2 = (df +  2.0D+00  * real ( icent, kind = dp )) /  2.0D+00
    lfact = gamma_log ( 1.0D+00 + dfd2 )
    lcntaj = dfd2 * log ( chid2 ) - chid2 - lfact
    centaj = exp ( lcntaj )
    sum1 = centwt * pcent
    !
    !  Sum backwards from the central term towards zero.
    !  Quit whenever either
    !  (1) the zero term is reached, or
    !  (2) the term gets small relative to the sum, or
    !  (3) More than NTIRED terms are totaled.
    !
    iterb = 0
    sumadj = 0.0D+00
    adj = centaj
    wt = centwt
    i = icent
    term = 0.0D+00

    do

        !dfd2 = dg(i) /  2.0D+00
        dfd2 = (df +  2.0D+00  * real ( i, kind = dp )) /  2.0D+00
        !
        !  Adjust chi-square for two fewer degrees of freedom.
        !  The adjusted value ends up in PTERM.
        !
        adj = adj * dfd2 / chid2
        sumadj = sumadj + adj
        pterm = pcent + sumadj
        !
        !  Adjust Poisson weight for J decreased by one.
        !
        wt = wt * ( i / xnonc )
        term = wt * pterm
        sum1 = sum1 + term
        i = i - 1
        iterb = iterb + 1

        !if ( ntired < iterb .or. qsmall ( term ) .or. i == 0 ) then
        if ( ntired < iterb .or. ( sum1 < 1.0D-20 .or. term < eps * sum1 ) .or. i == 0 ) then
            exit
        end if

    end do

    iterf = 0
    !
    !  Now sum forward from the central term towards infinity.
    !  Quit when either
    !    (1) the term gets small relative to the sum, or
    !    (2) More than NTIRED terms are totaled.
    !
    sumadj = centaj
    adj = centaj
    wt = centwt
    i = icent
    !
    !  Update weights for next higher J.
    !
    do

        wt = wt * ( xnonc / ( i + 1 ) )
        !
        !  Calculate PTERM and add term to sum.
        !
        pterm = pcent - sumadj
        term = wt * pterm
        sum1 = sum1 + term
        !
        !  Update adjustment term for DF for next iteration.
        !
        i = i + 1
        !dfd2 = dg(i) /  2.0D+00
        dfd2 = (df +  2.0D+00  * real ( i, kind = dp )) /  2.0D+00
        adj = adj * chid2 / dfd2
        sumadj = sumadj + adj
        iterf = iterf + 1

        !if ( ntired < iterf .or. qsmall ( term ) ) then
        if ( ntired < iterf .or. (sum1 < 1.0D-20 .or. term < eps * sum1) ) then
            exit
        end if

    end do

    cum = sum1
    ccum = 0.5D+00 + ( 0.5D+00 - cum )

    return
    end subroutine cumchn

    subroutine cumf ( f, dfn, dfd, cum, ccum )

    !*****************************************************************************80
    !
    !! CUMF evaluates the cumulative F distribution.
    !
    !  Discussion:
    !
    !    This routine computes the integral from 0 to F of the F density with DFN
    !    numerator and DFD denominator degrees of freedom.
    !
    !  Author:
    !
    !    Barry Brown, James Lovato, Kathy Russell
    !
    !  Reference:
    !
    !    Milton Abramowitz, Irene Stegun,
    !    Handbook of Mathematical Functions
    !    1966, Formula 26.5.28.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) F, the upper limit of integration.
    !
    !    Input, real ( kind = dp ) DFN, DFD, the number of degrees of
    !    freedom for the numerator and denominator.
    !
    !    Output, real ( kind = dp ) CUM, CCUM, the value of the F CDF and
    !    the complementary F CDF.
    !
    implicit none

    real ( kind = dp ) :: f,dfn,dfd,cum,ccum

    real ( kind = dp ) :: dsum,prod,xx, yy
    integer ( kind = int32 ) :: ierr

    if ( f <= 0.0D+00 ) then
        cum = 0.0D+00
        ccum = 1.0D+00
        return
    elseif (f /= f) then
        ! added by Patrick Macnamara to protect against
        ! scenario in which f is NaN
        cum = 0.0d0/0.0d0
        ccum = 0.0d0/0.0d0
        return
    end if

    prod = dfn * f
    !
    !  XX is such that the incomplete beta with parameters
    !  DFD/2 and DFN/2 evaluated at XX is 1 - CUM or CCUM
    !
    !  YY is 1 - XX
    !
    !  Calculate the smaller of XX and YY accurately.
    !
    dsum = dfd + prod
    xx = dfd / dsum

    if ( 0.5D+00 < xx ) then
        yy = prod / dsum
        xx = 1.0D+00 - yy
    else
        yy = 1.0D+00 - xx
    end if

    call beta_inc ( 0.5D+00 * dfd, 0.5D+00 * dfn, xx, yy, ccum, cum, ierr )

    return
    end subroutine cumf

    subroutine cumfnc ( f, dfn, dfd, pnonc, cum, ccum )

    !*****************************************************************************80
    !
    !! CUMFNC evaluates the cumulative noncentral F distribution.
    !
    !  Discussion:
    !
    !    This routine computes the noncentral F distribution with DFN and DFD
    !    degrees of freedom and noncentrality parameter PNONC.
    !
    !    The series is calculated backward and forward from J = LAMBDA/2
    !    (this is the term with the largest Poisson weight) until
    !    the convergence criterion is met.
    !
    !    The sum continues until a succeeding term is less than EPS
    !    times the sum or the sum is very small.  EPS is
    !    set to 1.0D-4 in a data statement which can be changed.
    !
    !    The original version of this routine allowed the input values
    !    of DFN and DFD to be negative (nonsensical) or zero (which
    !    caused numerical overflow.)  I have forced both these values
    !    to be at least 1.
    !
    !  Modified:
    !
    !    19 May 2007
    !
    !  Author:
    !
    !    Barry Brown, James Lovato, Kathy Russell
    !
    !  Reference:
    !
    !    Milton Abramowitz, Irene Stegun,
    !    Handbook of Mathematical Functions
    !    1966, Formula 26.5.16, 26.6.17, 26.6.18, 26.6.20.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) F, the upper limit of integration.
    !
    !    Input, real ( kind = dp ) DFN, DFD, the number of degrees of freedom
    !    in the numerator and denominator.  Both DFN and DFD must be positive,
    !    and normally would be integers.  This routine requires that they
    !    be no less than 1.
    !
    !    Input, real ( kind = dp ) PNONC, the noncentrality parameter.
    !
    !    Output, real ( kind = dp ) CUM, CCUM, the noncentral F CDF and
    !    complementary CDF.
    !
    implicit none

    real ( kind = dp ) :: f,dfn,dfd,pnonc,cum,ccum

    real ( kind = dp ) :: adn,arg1,aup,b,betdn,betup,centwt,dnterm,dsum,dummy
    real ( kind = dp ) :: expon,prod,sum1,upterm,xmult,xnonc,xx,yy
    integer ( kind = int32 ) :: i,icent,ierr

    real ( kind = dp ), parameter :: eps = 0.0001D+00

    if ( f <= 0.0D+00 ) then
        cum = 0.0D+00
        ccum = 1.0D+00
        return
    end if

    if ( dfn < 1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CUMFNC - Fatal error!'
        write ( *, '(a)' ) '  DFN < 1.'
        stop
    end if

    if ( dfd < 1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CUMFNC - Fatal error!'
        write ( *, '(a)' ) '  DFD < 1.'
        stop
    end if
    !
    !  Handle case in which the noncentrality parameter is essentially zero.
    !
    if ( pnonc < 1.0D-10 ) then
        call cumf ( f, dfn, dfd, cum, ccum )
        return
    end if

    xnonc = pnonc /  2.0D+00
    !
    !  Calculate the central term of the Poisson weighting factor.
    !
    icent = int ( xnonc )

    if ( icent == 0 ) then
        icent = 1
    end if
    !
    !  Compute central weight term.
    !
    centwt = exp ( -xnonc + icent * log ( xnonc ) &
        - gamma_log ( real ( icent + 1, kind = dp  ) ) )
    !
    !  Compute central incomplete beta term.
    !  Ensure that minimum of arg to beta and 1 - arg is computed accurately.
    !
    prod = dfn * f
    dsum = dfd + prod
    yy = dfd / dsum

    if ( 0.5D+00 < yy ) then
        xx = prod / dsum
        yy = 1.0D+00 - xx
    else
        xx = 1.0D+00 - yy
    end if

    arg1 = 0.5D+00 * dfn + real ( icent, kind = dp )
    call beta_inc ( arg1, 0.5D+00*dfd, xx, yy, betdn, dummy, ierr )

    adn = dfn / 2.0D+00 + real ( icent, kind = dp )
    aup = adn
    b = dfd / 2.0D+00
    betup = betdn
    sum1 = centwt * betdn
    !
    !  Now sum terms backward from ICENT until convergence or all done.
    !
    xmult = centwt
    i = icent
    dnterm = exp ( gamma_log ( adn + b ) &
        - gamma_log ( adn + 1.0D+00 ) &
        - gamma_log ( b ) + adn * log ( xx ) + b * log ( yy ) )

    do

        if ( i <= 0 ) then
            exit
        end if

        if ( sum1 < epsilon ( xmult * betdn ) .or. &
            xmult * betdn < eps * sum1 ) then
        exit
        end if

        xmult = xmult * ( real ( i, kind = dp ) / xnonc )
        i = i - 1
        adn = adn - 1.0D+00
        dnterm = ( adn + 1.0D+00 ) / ( ( adn + b ) * xx ) * dnterm
        betdn = betdn + dnterm
        sum1 = sum1 + xmult * betdn

    end do

    i = icent + 1
    !
    !  Now sum forward until convergence.
    !
    xmult = centwt

    if ( ( aup - 1.0D+00 + b ) == 0 ) then

        expon = - gamma_log ( aup ) - gamma_log ( b ) &
            + ( aup - 1.0D+00 ) * log ( xx ) + b * log ( yy )

    else

        expon = gamma_log ( aup - 1.0D+00 + b ) - gamma_log ( aup ) &
            - gamma_log ( b ) + ( aup - 1.0D+00 ) * log ( xx ) + b * log ( yy )

    end if
    !
    !  The fact that DCDFLIB blithely assumes that 1.0E+30 is a reasonable
    !  value to plug into any function, and that G95 computes corresponding
    !  function values of, say 1.0E-303, and then chokes with a floating point
    !  error when asked to combine such a value with a reasonable floating
    !  point quantity, has driven me to the following sort of check that
    !  was last fashionable in the 1960's!
    !
    if ( expon <= log ( epsilon ( expon ) ) ) then
        upterm = 0.0D+00
    else
        upterm = exp ( expon )
    end if

    do

        xmult = xmult * ( xnonc / real ( i, kind = dp ) )
        i = i + 1
        aup = aup + 1.0D+00
        upterm = ( aup + b -  2.0D+00  ) * xx / ( aup - 1.0D+00 ) * upterm
        betup = betup - upterm
        sum1 = sum1 + xmult * betup

        if ( sum1 < epsilon ( xmult * betup ) .or. xmult * betup < eps * sum1 ) then
            exit
        end if

    end do

    cum = sum1
    ccum = 0.5D+00 + ( 0.5D+00 - cum )

    return
    end subroutine cumfnc

    subroutine cumgam ( x, a, cum, ccum )

    !*****************************************************************************80
    !
    !! CUMGAM evaluates the cumulative incomplete gamma distribution.
    !
    !  Discussion:
    !
    !    This routine computes the cumulative distribution function of the
    !    incomplete gamma distribution, i.e., the integral from 0 to X of
    !
    !      (1/GAM(A))*EXP(-T)*T^(A-1) DT
    !
    !    where GAM(A) is the complete gamma function of A:
    !
    !      GAM(A) = integral from 0 to infinity of EXP(-T)*T^(A-1) DT
    !
    !  Author:
    !
    !    Barry Brown, James Lovato, Kathy Russell
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) X, the upper limit of integration.
    !
    !    Input, real ( kind = dp ) A, the shape parameter of the incomplete
    !    Gamma distribution.
    !
    !    Output, real ( kind = dp ) CUM, CCUM, the incomplete Gamma CDF and
    !    complementary CDF.
    !
    implicit none

    real ( kind = dp ) :: x,a,cum,ccum

    if ( x <= 0.0D+00 ) then

        cum = 0.0D+00
        ccum = 1.0D+00

    else

        call gamma_inc ( a, x, cum, ccum, 0 )

    end if

    return
    end subroutine cumgam

    subroutine cumnbn ( f, s, pr, ompr, cum, ccum )

    !*****************************************************************************80
    !
    !! CUMNBN evaluates the cumulative negative binomial distribution.
    !
    !  Discussion:
    !
    !    This routine returns the probability that there will be F or
    !    fewer failures before there are S successes, with each binomial
    !    trial having a probability of success PR.
    !
    !    Prob(# failures = F | S successes, PR)  =
    !                        ( S + F - 1 )
    !                        (            ) * PR^S * (1-PR)^F
    !                        (      F     )
    !
    !  Author:
    !
    !    Barry Brown, James Lovato, Kathy Russell
    !
    !  Reference:
    !
    !    Milton Abramowitz, Irene Stegun,
    !    Handbook of Mathematical Functions
    !    1966, Formula 26.5.26.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) F, the number of failures.
    !
    !    Input, real ( kind = dp ) S, the number of successes.
    !
    !    Input, real ( kind = dp ) PR, OMPR, the probability of success on
    !    each binomial trial, and the value of (1-PR).
    !
    !    Output, real ( kind = dp ) CUM, CCUM, the negative binomial CDF,
    !    and the complementary CDF.
    !
    implicit none

    real ( kind = dp ) :: f,s,pr,ompr,cum,ccum

    call cumbet ( pr, ompr, s, f+1.D+00, cum, ccum )

    return
    end subroutine cumnbn

    function normcdf ( x, upper )

    !*****************************************************************************80
    !
    !! NORMCDF computes the cumulative density of the standard normal distribution.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) x, is one endpoint of the semi-infinite interval
    !    over which the integration takes place.
    !
    !    Input, logical UPPER, determines whether the upper or lower
    !    interval is to be integrated:
    !    .TRUE.  => integrate from X to + Infinity (i.e., complementary CDF)
    !    .FALSE. => integrate from - Infinity to X (i.e., CDF)
    !
    !    Output, real ( kind = dp ) NORMCDF, the integral of the standard normal
    !    distribution over the desired interval.
    !
    implicit none

    real ( kind = dp ) :: x
    logical :: upper
    real ( kind = dp ) :: normcdf

    ! locals
    real( kind = dp ) :: cum, ccum

    ! evaluate normal cdf
    call cumnor(x, cum, ccum)

    ! determine which output to return
    normcdf = merge(ccum, cum, upper)

    end function normcdf

    subroutine cumnor ( arg, cum, ccum )

    !*****************************************************************************80
    !
    !! CUMNOR computes the cumulative normal distribution.
    !
    !  Discussion:
    !
    !    This function evaluates the normal distribution function:
    !
    !                              / x
    !                     1       |       -t*t/2
    !          P(x) = ----------- |      e       dt
    !                 sqrt(2 pi)  |
    !                             /-oo
    !
    !    This transportable program uses rational functions that
    !    theoretically approximate the normal distribution function to
    !    at least 18 significant decimal digits.  The accuracy achieved
    !    depends on the arithmetic system, the compiler, the intrinsic
    !    functions, and proper selection of the machine dependent
    !    constants.
    !
    !  Author:
    !
    !    William Cody
    !    Mathematics and Computer Science Division
    !    Argonne National Laboratory
    !    Argonne, IL 60439
    !
    !  Reference:
    !
    !    William Cody,
    !    Rational Chebyshev approximations for the error function,
    !    Mathematics of Computation,
    !    1969, pages 631-637.
    !
    !    William Cody,
    !    Algorithm 715:
    !    SPECFUN - A Portable FORTRAN Package of Special Function Routines
    !    and Test Drivers,
    !    ACM Transactions on Mathematical Software,
    !    Volume 19, Number 1, 1993, pages 22-32.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) ARG, the upper limit of integration.
    !
    !    Output, real ( kind = dp ) CUM, CCUM, the Normal density CDF and
    !    complementary CDF.
    !
    !  Local Parameters:
    !
    !    Local, real ( kind = dp ) EPS, the argument below which anorm(x)
    !    may be represented by 0.5 and above which  x*x  will not underflow.
    !    A conservative value is the largest machine number X
    !    such that   1.0D+00 + X = 1.0D+00   to machine precision.
    !
    implicit none

    real ( kind = dp ) :: arg,cum,ccum

    real ( kind = dp ) :: del,eps,temp,x,xden,xnum,y,xsq
    integer ( kind = int32 ) :: i

    real ( kind = dp ), parameter, dimension ( 5 ) :: a = (/ &
        2.2352520354606839287D+00, &
        1.6102823106855587881D+02, &
        1.0676894854603709582D+03, &
        1.8154981253343561249D+04, &
        6.5682337918207449113D-02 /)
    real ( kind = dp ), parameter, dimension ( 4 ) :: b = (/ &
        4.7202581904688241870D+01, &
        9.7609855173777669322D+02, &
        1.0260932208618978205D+04, &
        4.5507789335026729956D+04 /)
    real ( kind = dp ), parameter, dimension ( 9 ) :: c = (/ &
        3.9894151208813466764D-01, &
        8.8831497943883759412D+00, &
        9.3506656132177855979D+01, &
        5.9727027639480026226D+02, &
        2.4945375852903726711D+03, &
        6.8481904505362823326D+03, &
        1.1602651437647350124D+04, &
        9.8427148383839780218D+03, &
        1.0765576773720192317D-08 /)
    real ( kind = dp ), parameter, dimension ( 8 ) :: d = (/ &
        2.2266688044328115691D+01, &
        2.3538790178262499861D+02, &
        1.5193775994075548050D+03, &
        6.4855582982667607550D+03, &
        1.8615571640885098091D+04, &
        3.4900952721145977266D+04, &
        3.8912003286093271411D+04, &
        1.9685429676859990727D+04 /)
    real ( kind = dp ), parameter, dimension ( 6 ) :: p = (/ &
        2.1589853405795699D-01, &
        1.274011611602473639D-01, &
        2.2235277870649807D-02, &
        1.421619193227893466D-03, &
        2.9112874951168792D-05, &
        2.307344176494017303D-02 /)
    real ( kind = dp ), parameter, dimension ( 5 ) :: q = (/ &
        1.28426009614491121D+00, &
        4.68238212480865118D-01, &
        6.59881378689285515D-02, &
        3.78239633202758244D-03, &
        7.29751555083966205D-05 /)

    real ( kind = dp ), parameter :: root32 = 5.656854248D+00
    real ( kind = dp ), parameter :: sixten = 16.0D+00
    real ( kind = dp ), parameter :: sqrpi = 3.9894228040143267794D-01
    real ( kind = dp ), parameter :: thrsh = 0.66291D+00

    !
    !  Machine dependent constants
    !
    eps = epsilon ( 1.0D+00 ) * 0.5D+00

    x = arg
    y = abs ( x )

    if ( y <= thrsh ) then
        !
        !  Evaluate  anorm  for  |X| <= 0.66291
        !
        if ( eps < y ) then
            xsq = x * x
        else
            xsq = 0.0D+00
        end if

        xnum = a(5) * xsq
        xden = xsq
        do i = 1, 3
            xnum = ( xnum + a(i) ) * xsq
            xden = ( xden + b(i) ) * xsq
        end do
        cum = x * ( xnum + a(4) ) / ( xden + b(4) )
        temp = cum
        cum = 0.5D+00 + temp
        ccum = 0.5D+00 - temp
        !
        !  Evaluate ANORM for 0.66291 <= |X| <= sqrt(32)
        !
    else if ( y <= root32 ) then

        xnum = c(9) * y
        xden = y
        do i = 1, 7
            xnum = ( xnum + c(i) ) * y
            xden = ( xden + d(i) ) * y
        end do
        cum = ( xnum + c(8) ) / ( xden + d(8) )
        xsq = aint ( y * sixten ) / sixten
        del = ( y - xsq ) * ( y + xsq )
        cum = exp ( - xsq * xsq * 0.5D+00 ) * exp ( -del * 0.5D+00 ) * cum
        ccum = 1.0D+00 - cum

        if ( 0.0D+00 < x ) then
            call r8_swap ( cum, ccum )
        end if
        !
        !  Evaluate ANORM for sqrt(32) < |X|.
        !
    else

        cum = 0.0D+00
        xsq = 1.0D+00 / ( x * x )
        xnum = p(6) * xsq
        xden = xsq
        do i = 1, 4
            xnum = ( xnum + p(i) ) * xsq
            xden = ( xden + q(i) ) * xsq
        end do

        cum = xsq * ( xnum + p(5) ) / ( xden + q(5) )
        cum = ( sqrpi - cum ) / y
        xsq = aint ( x * sixten ) / sixten
        del = ( x - xsq ) * ( x + xsq )
        cum = exp ( - xsq * xsq * 0.5D+00 ) &
            * exp ( - del * 0.5D+00 ) * cum
        ccum = 1.0D+00 - cum

        if ( 0.0D+00 < x ) then
            call r8_swap ( cum, ccum )
        end if

    end if

    if ( cum < tiny ( cum ) ) then
        cum = 0.0D+00
    end if

    if ( ccum < tiny ( ccum ) ) then
        ccum = 0.0D+00
    end if

    return
    end subroutine cumnor

    subroutine cumpoi ( s, xlam, cum, ccum )

    !*****************************************************************************80
    !
    !! CUMPOI evaluates the cumulative Poisson distribution.
    !
    !  Discussion:
    !
    !    This routine returns the probability of S or fewer events in a Poisson
    !    distribution with mean XLAM.
    !
    !  Author:
    !
    !    Barry Brown, James Lovato, Kathy Russell
    !
    !  Reference:
    !
    !    Milton Abramowitz, Irene Stegun,
    !    Handbook of Mathematical Functions,
    !    Formula 26.4.21.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) S, the upper limit of cumulation of the
    !    Poisson density function.
    !
    !    Input, real ( kind = dp ) XLAM, the mean of the Poisson distribution.
    !
    !    Output, real ( kind = dp ) CUM, CCUM, the Poisson density CDF and
    !    complementary CDF.
    !
    implicit none

    real ( kind = dp ) :: s,xlam,cum,ccum

    real ( kind = dp ) :: chi,df

    df =  2.0D+00  * ( s + 1.0D+00 )
    chi =  2.0D+00  * xlam

    call cumchi ( chi, df, ccum, cum )

    return
    end subroutine cumpoi

    subroutine cumt ( t, df, cum, ccum )

    !*****************************************************************************80
    !
    !! CUMT evaluates the cumulative T distribution.
    !
    !  Author:
    !
    !    Barry Brown, James Lovato, Kathy Russell
    !
    !  Reference:
    !
    !    Milton Abramowitz, Irene Stegun,
    !    Handbook of Mathematical Functions,
    !    Formula 26.5.27.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) T, the upper limit of integration.
    !
    !    Input, real ( kind = dp ) DF, the number of degrees of freedom of
    !    the T distribution.
    !
    !    Output, real ( kind = dp ) CUM, CCUM, the T distribution CDF and
    !    complementary CDF.
    !
    implicit none

    real ( kind = dp ) :: t,df,cum,ccum

    real ( kind = dp ) :: a,oma,xx,yy

    xx = df / ( df + t**2 )
    yy = t**2 / ( df + t**2 )

    call cumbet ( xx, yy, 0.5D+00*df, 0.5D+00, a, oma )

    if ( t /= t ) then
        ! added by Patrick Macnamara to protect
        ! against scenario in which t is NaN
        cum = 0.0d0/0.0d0
        ccum = 0.0d0/0.0d0
    elseif ( t <= 0.0D+00 ) then
        cum = 0.5D+00 * a
        ccum = oma + cum
    else
        ccum = 0.5D+00 * a
        cum = oma + ccum
    end if

    return
    end subroutine cumt

    function dbetrm ( a, b )

    !*****************************************************************************80
    !
    !! DBETRM computes the Sterling remainder for the complete beta function.
    !
    !  Discussion:
    !
    !    Log(Beta(A,B)) = Lgamma(A) + Lgamma(B) - Lgamma(A+B)
    !    where Lgamma is the log of the (complete) gamma function
    !
    !    Let ZZ be approximation obtained if each log gamma is approximated
    !    by Sterling's formula, i.e.,
    !
    !      Sterling(Z) = log ( sqrt ( 2 * PI ) ) + ( Z - 0.5 ) * log ( Z ) - Z
    !
    !    The Sterling remainder is Log(Beta(A,B)) - ZZ.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) A, B, the parameters of the Beta function.
    !
    !    Output, real ( kind = dp ) DBETRM, the Sterling remainder.
    !
    implicit none

    real ( kind = dp ) :: a,b
    real ( kind = dp ) :: dbetrm

    !
    !  Try to sum from smallest to largest.
    !
    dbetrm = -dstrem ( a + b )
    dbetrm = dbetrm + dstrem ( max ( a, b ) )
    dbetrm = dbetrm + dstrem ( min ( a, b ) )

    return
    end function dbetrm

    function dexpm1 ( x )

    !*****************************************************************************80
    !
    !! DEXPM1 evaluates the function EXP(X) - 1.
    !
    !  Author:
    !
    !    Armido DiDinato, Alfred Morris
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) X, the value at which exp(X)-1 is desired.
    !
    !    Output, real ( kind = dp ) DEXPM1, the value of exp(X)-1.
    !
    implicit none

    real ( kind = dp ) :: x
    real ( kind = dp ) :: dexpm1

    real ( kind = dp ) :: bot,top,w
    real ( kind = dp ), parameter :: p1 =  0.914041914819518D-09
    real ( kind = dp ), parameter :: p2 =  0.238082361044469D-01
    real ( kind = dp ), parameter :: q1 = -0.499999999085958D+00
    real ( kind = dp ), parameter :: q2 =  0.107141568980644D+00
    real ( kind = dp ), parameter :: q3 = -0.119041179760821D-01
    real ( kind = dp ), parameter :: q4 =  0.595130811860248D-03

    if ( abs ( x ) <= 0.15D+00 ) then

        top = ( p2 * x + p1 ) * x + 1.0D+00
        bot = ((( q4 * x + q3 ) * x + q2 ) * x + q1 ) * x + 1.0D+00
        dexpm1 = x * ( top / bot )

    else

        w = exp ( x )

        if ( x <= 0.0D+00 ) then
            dexpm1 = ( w - 0.5D+00 ) - 0.5D+00
        else
            dexpm1 = w * ( 0.5D+00 &
                + ( 0.5D+00 - 1.0D+00 / w ))
        end if

    end if

    return
    end function dexpm1

    function norminv ( p )
    !*****************************************************************************80
    !
    !! NORMINV computes the inverse of the normal distribution.
    !
    !  Discussion:
    !
    !    This routine returns X such that
    !
    !      NORMCDF(X) = P,
    !
    !    that is, so that
    !
    !      P = integral ( -oo <= T <= X ) exp(-U*U/2)/sqrt(2*PI) dU
    !
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) P, the probability
    !
    !    Output, real ( kind = dp ) NORMINV, the argument X for which the
    !    Normal CDF has the value P.
    !
    implicit none

    real ( kind = dp ) :: p
    real ( kind = dp ) :: norminv

    ! locals
    real ( kind = dp ) :: q

    ! determine q
    q = 1.0d0 - p

    norminv = dinvnr(p,q)

    end function norminv

    function dinvnr ( p, q )

    !*****************************************************************************80
    !
    !! DINVNR computes the inverse of the normal distribution.
    !
    !  Discussion:
    !
    !    This routine returns X such that
    !
    !      CUMNOR(X) = P,
    !
    !    that is, so that
    !
    !      P = integral ( -oo <= T <= X ) exp(-U*U/2)/sqrt(2*PI) dU
    !
    !    The rational function on page 95 of Kennedy and Gentle is used as a
    !    starting value for the Newton method of finding roots.
    !
    !  Reference:
    !
    !    William Kennedy, James Gentle,
    !    Statistical Computing,
    !    Marcel Dekker, NY, 1980,
    !    QA276.4 K46
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) P, Q, the probability, and the complementary
    !    probability.
    !
    !    Output, real ( kind = dp ) DINVNR, the argument X for which the
    !    Normal CDF has the value P.
    !
    implicit none

    real ( kind = dp ) :: p,q
    real ( kind = dp ) :: dinvnr

    real ( kind = dp ) :: ccum,cum,dx,pp,strtx,xcur
    integer ( kind = int32 ) :: i

    real ( kind = dp ), parameter :: eps = 1.0D-13
    real ( kind = dp ), parameter :: r2pi = 0.3989422804014326D+00
    integer ( kind = int32 ), parameter :: maxit = 100

    pp = min ( p, q )
    strtx = stvaln ( pp )
    xcur = strtx
    !
    !  Newton iterations.
    !
    do i = 1, maxit

        call cumnor ( xcur, cum, ccum )
        dx = ( cum - pp ) / ( r2pi * exp ( -0.5D+00 * xcur * xcur ) )
        xcur = xcur - dx

        if ( abs ( dx / xcur ) < eps ) then
            if ( p <= q ) then
                dinvnr = xcur
            else
                dinvnr = -xcur
            end if
            return
        end if

    end do

    if ( p <= q ) then
        dinvnr = strtx
    else
        dinvnr = -strtx
    end if

    return
    end function dinvnr

    function dlanor ( x )

    !*****************************************************************************80
    !
    !! DLANOR evaluates the logarithm of the asymptotic Normal CDF.
    !
    !  Discussion:
    !
    !    This routine computes the logarithm of the cumulative normal distribution
    !    from abs ( x ) to infinity for  5 <= abs ( X ).
    !
    !    The relative error at X = 5 is about 0.5D-5.
    !
    !  Reference:
    !
    !    Milton Abramowitz, Irene Stegun,
    !    Handbook of Mathematical Functions
    !    1966, Formula 26.2.12.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) X, the value at which the Normal CDF is to be
    !    evaluated.  It is assumed that 5 <= abs ( X ).
    !
    !    Output, real ( kind = dp ) DLANOR, the logarithm of the asymptotic
    !    Normal CDF.
    !
    implicit none

    real ( kind = dp ) :: x
    real ( kind = dp ) :: dlanor

    real ( kind = dp ) :: approx, correc,xx,xx2
    real ( kind = dp ), save, dimension ( 0:11 ) :: coef = (/ &
        -1.0D+00,  3.0D+00,  -15.0D+00,  105.0D+00,  -945.0D+00,  &
        10395.0D+00, -135135.0D+00,  2027025.0D+00,  -34459425.0D+00, &
        654729075.0D+00, -13749310575D+00,  316234143225.0D+00 /)
    real ( kind = dp ), parameter :: dlsqpi = 0.91893853320467274177D+00

    xx = abs ( x )

    if ( abs ( x ) < 5.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DLANOR - Fatal error!'
        write ( *, '(a)' ) '  The argument X is too small.'
    end if

    approx = - dlsqpi - 0.5D+00 * x**2 - log ( abs ( x ) )

    xx2 = xx * xx
    correc = eval_pol ( coef, 11, 1.0D+00 / xx2 ) / xx2
    correc = alnrel ( correc )

    dlanor = approx + correc

    return
    end function dlanor

    function dstrem ( z )

    !*****************************************************************************80
    !
    !! DSTREM computes the Sterling remainder ln ( Gamma ( Z ) ) - Sterling ( Z ).
    !
    !  Discussion:
    !
    !    This routine returns
    !
    !      ln ( Gamma ( Z ) ) - Sterling ( Z )
    !
    !    where Sterling(Z) is Sterling's approximation to ln ( Gamma ( Z ) ).
    !
    !    Sterling(Z) = ln ( sqrt ( 2 * PI ) ) + ( Z - 0.5 ) * ln ( Z ) - Z
    !
    !    If 6 <= Z, the routine uses 9 terms of a series in Bernoulli numbers,
    !    with values calculated using Maple.
    !
    !    Otherwise, the difference is computed explicitly.
    !
    !  Modified:
    !
    !    14 June 2004
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) Z, the value at which the Sterling
    !    remainder is to be calculated.  Z must be positive.
    !
    !    Output, real ( kind = dp ) DSTREM, the Sterling remainder.
    !
    implicit none

    real ( kind = dp ) :: z
    real ( kind = dp ) :: dstrem

    real ( kind = dp ), parameter :: hln2pi = 0.91893853320467274178D+00
    integer ( kind = int32 ), parameter :: ncoef = 9

    real ( kind = dp ) :: sterl
    real ( kind = dp ), parameter, dimension ( 0:ncoef ) :: coef = (/ &
        0.0D+00, &
        0.0833333333333333333333333333333D+00, &
        -0.00277777777777777777777777777778D+00, &
        0.000793650793650793650793650793651D+00, &
        -0.000595238095238095238095238095238D+00, &
        0.000841750841750841750841750841751D+00, &
        -0.00191752691752691752691752691753D+00, &
        0.00641025641025641025641025641026D+00, &
        -0.0295506535947712418300653594771D+00, &
        0.179644372368830573164938490016D+00 /)

    if ( z <= 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DSTREM - Fatal error!'
        write ( *, '(a)' ) '  Zero or negative argument Z.'
        stop
    end if

    if ( 6.0D+00 < z ) then
        dstrem = eval_pol ( coef, ncoef, 1.0D+00 / z**2 ) * z
    else
        sterl = hln2pi + ( z - 0.5D+00 ) * log ( z ) - z
        dstrem = gamma_log ( z ) - sterl
    end if

    return
    end function dstrem

    function dt1 ( p, q, df )

    !*****************************************************************************80
    !
    !! DT1 computes an approximate inverse of the cumulative T distribution.
    !
    !  Discussion:
    !
    !    This routine returns the inverse of the T distribution function, that is,
    !    the integral from 0 to INVT of the T density is P.  This is an
    !    initial approximation.
    !
    !    Thanks to Charles Katholi for pointing out that the RESHAPE
    !    function should not use a range in the "SHAPE" field (0:4,4),
    !    but simply the number of rows and columns (5,4), JVB, 04 May 2006.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) P, Q, the value whose inverse from the
    !    T distribution CDF is desired, and the value (1-P).
    !
    !    Input, real ( kind = dp ) DF, the number of degrees of freedom of the
    !    T distribution.
    !
    !    Output, real ( kind = dp ) DT1, the approximate value of X for which
    !    the T density CDF with DF degrees of freedom has value P.
    !
    implicit none

    real ( kind = dp ) :: p,q,df
    real ( kind = dp ) :: dt1

    real ( kind = dp ) :: denpow,sum1,term,x,xp,xx
    integer ( kind = int32 ) :: i

    real ( kind = dp ), dimension(0:4,4) :: coef = reshape ( (/ &
        1.0D+00,     1.0D+00,    0.0D+00,   0.0D+00,  0.0D+00, &
        3.0D+00,    16.0D+00,    5.0D+00,   0.0D+00,  0.0D+00, &
        -15.0D+00,    17.0D+00,   19.0D+00,   3.0D+00,  0.0D+00, &
        -945.0D+00, -1920.0D+00, 1482.0D+00, 776.0D+00, 79.0D+00/), (/ 5, 4 /) )
    real ( kind = dp ), parameter, dimension ( 4 ) :: denom = (/ &
        4.0D+00, 96.0D+00, 384.0D+00, 92160.0D+00 /)
    integer ( kind = int32 ), parameter, dimension ( 4 ) :: ideg = (/ 1, 2, 3, 4 /)

    x = abs ( dinvnr ( p, q ) )
    xx = x * x

    sum1 = x
    denpow = 1.0D+00
    do i = 1, 4
        term = eval_pol ( coef(0,i), ideg(i), xx ) * x
        denpow = denpow * df
        sum1 = sum1 + term / ( denpow * denom(i) )
    end do

    if ( 0.5D+00 <= p ) then
        xp = sum1
    else
        xp = -sum1
    end if

    dt1 = xp

    return
    end function dt1

    subroutine erf_values ( n_data, x, fx )

    !*****************************************************************************80
    !
    !! ERF_VALUES returns some values of the ERF or "error" function.
    !
    !  Discussion:
    !
    !    ERF(X) = ( 2 / sqrt ( PI ) * integral ( 0 <= T <= X ) exp ( - T^2 ) dT
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    17 April 2001
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Milton Abramowitz, Irene Stegun,
    !    Handbook of Mathematical Functions,
    !    US Department of Commerce, 1964.
    !
    !  Parameters:
    !
    !    Input/output, integer ( kind = int32 ) N_DATA.  The user sets N_DATA to 0
    !    before the first call.  On each call, the routine increments N_DATA by 1,
    !    and returns the corresponding data; when there is no more data, the
    !    output value of N_DATA will be 0 again.
    !
    !    Output, real ( kind = dp ) X, the argument of the function.
    !
    !    Output, real ( kind = dp ) FX, the value of the function.
    !
    implicit none

    integer ( kind = int32 ) :: n_data
    real ( kind = dp ) :: x,fx

    integer ( kind = int32 ), parameter :: n_max = 21
    real ( kind = dp ), save, dimension ( n_max ) :: fx_vec = (/ &
        0.0000000000D+00, 0.1124629160D+00, 0.2227025892D+00, 0.3286267595D+00, &
        0.4283923550D+00, 0.5204998778D+00, 0.6038560908D+00, 0.6778011938D+00, &
        0.7421009647D+00, 0.7969082124D+00, 0.8427007929D+00, 0.8802050696D+00, &
        0.9103139782D+00, 0.9340079449D+00, 0.9522851198D+00, 0.9661051465D+00, &
        0.9763483833D+00, 0.9837904586D+00, 0.9890905016D+00, 0.9927904292D+00, &
        0.9953222650D+00 /)
    real ( kind = dp ), save, dimension ( n_max ) :: x_vec = (/ &
        0.0D+00, 0.1D+00, 0.2D+00, 0.3D+00, &
        0.4D+00, 0.5D+00, 0.6D+00, 0.7D+00, &
        0.8D+00, 0.9D+00, 1.0D+00, 1.1D+00, &
        1.2D+00, 1.3D+00, 1.4D+00, 1.5D+00, &
        1.6D+00, 1.7D+00, 1.8D+00, 1.9D+00, &
        2.0D+00 /)

    if ( n_data < 0 ) then
        n_data = 0
    end if

    n_data = n_data + 1

    if ( n_max < n_data ) then
        n_data = 0
        x = 0.0D+00
        fx = 0.0D+00
    else
        x = x_vec(n_data)
        fx = fx_vec(n_data)
    end if

    return
    end subroutine erf_values

    function error_f ( x )

    !*****************************************************************************80
    !
    !! ERROR_F evaluates the error function.
    !
    !  Discussion:
    !
    !    Since some compilers already supply a routine named ERF which evaluates
    !    the error function, this routine has been given a distinct, if
    !    somewhat unnatural, name.
    !
    !    The function is defined by:
    !
    !      ERF(X) = ( 2 / sqrt ( PI ) )
    !        * Integral ( 0 <= T <= X ) EXP ( - T**2 ) dT.
    !
    !    Properties of the function include:
    !
    !      Limit ( X -> -Infinity ) ERF(X) =          -1.0;
    !                               ERF(0) =           0.0;
    !                               ERF(0.476936...) = 0.5;
    !      Limit ( X -> +Infinity ) ERF(X) =          +1.0.
    !
    !      0.5D+00 * ( ERF(X/sqrt(2)) + 1 ) = Normal_01_CDF(X)
    !
    !  Modified:
    !
    !    17 November 2006
    !
    !  Author:
    !
    !    Armido DiDinato, Alfred Morris
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) X, the argument.
    !
    !    Output, real ( kind = dp ) ERF, the value of the error function at X.
    !
    implicit none

    real ( kind = dp ) :: x
    real ( kind = dp ) :: error_f

    real ( kind = dp ) :: ax,bot,t,top,x2
    real ( kind = dp ), parameter :: c = 0.564189583547756D+00

    real ( kind = dp ), parameter, dimension ( 5 ) :: a = (/ &
        0.771058495001320D-04, &
        -0.133733772997339D-02, &
        0.323076579225834D-01, &
        0.479137145607681D-01, &
        0.128379167095513D+00 /)
    real ( kind = dp ), parameter, dimension ( 3 ) :: b = (/ &
        0.301048631703895D-02, &
        0.538971687740286D-01, &
        0.375795757275549D+00 /)
    real ( kind = dp ), dimension ( 8 ) :: p = (/   &
        -1.36864857382717D-07, 5.64195517478974D-01, &
        7.21175825088309D+00, 4.31622272220567D+01, &
        1.52989285046940D+02, 3.39320816734344D+02, &
        4.51918953711873D+02, 3.00459261020162D+02 /)
    real ( kind = dp ), dimension ( 8 ) :: q = (/ &
        1.00000000000000D+00, 1.27827273196294D+01, &
        7.70001529352295D+01, 2.77585444743988D+02, &
        6.38980264465631D+02, 9.31354094850610D+02, &
        7.90950925327898D+02, 3.00459260956983D+02 /)
    real ( kind = dp ), dimension ( 5 ) :: r = (/ &
        2.10144126479064D+00, 2.62370141675169D+01, &
        2.13688200555087D+01, 4.65807828718470D+00, &
        2.82094791773523D-01 /)
    real ( kind = dp ), parameter, dimension ( 4 ) :: s = (/ &
        9.41537750555460D+01, 1.87114811799590D+02, &
        9.90191814623914D+01, 1.80124575948747D+02 /)

    ax = abs ( x )

    if ( ax <= 0.5D+00 ) then

        t = x * x

        top = (((( a(1)   * t &
            + a(2) ) * t &
            + a(3) ) * t &
            + a(4) ) * t &
            + a(5) ) + 1.0D+00

        bot = (( b(1) * t + b(2) ) * t + b(3) ) * t + 1.0D+00
        error_f = ax * ( top / bot )

    else if ( ax <= 4.0D+00 ) then

        top = (((((( p(1)   * ax &
            + p(2) ) * ax &
            + p(3) ) * ax &
            + p(4) ) * ax &
            + p(5) ) * ax &
            + p(6) ) * ax &
            + p(7) ) * ax &
            + p(8)

        bot = (((((( q(1) * ax + q(2) ) * ax + q(3) ) * ax + q(4) ) * ax &
            + q(5) ) * ax + q(6) ) * ax + q(7) ) * ax + q(8)

        error_f = 0.5D+00 &
            + ( 0.5D+00 - exp ( - x * x ) * top / bot )

    else if ( ax < 5.8D+00 ) then

        x2 = x * x
        t = 1.0D+00 / x2

        top = ((( r(1) * t + r(2) ) * t + r(3) ) * t + r(4) ) * t + r(5)

        bot = ((( s(1) * t + s(2) ) * t + s(3) ) * t + s(4) ) * t &
            + 1.0D+00

        error_f = ( c - top / ( x2 * bot )) / ax
        error_f = 0.5D+00 &
            + ( 0.5D+00 - exp ( - x2 ) * error_f )

    else

        error_f = 1.0D+00

    end if

    if ( x < 0.0D+00 ) then
        error_f = -error_f
    end if

    return
    end function error_f

    function error_fc ( ind, x )

    !*****************************************************************************80
    !
    !! ERROR_FC evaluates the complementary error function.
    !
    !  Modified:
    !
    !    09 December 1999
    !
    !  Author:
    !
    !    Armido DiDinato, Alfred Morris
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = int32 ) IND, chooses the scaling.
    !    If IND is nonzero, then the value returned has been multiplied by
    !    EXP(X*X).
    !
    !    Input, real ( kind = dp ) X, the argument of the function.
    !
    !    Output, real ( kind = dp ) ERROR_FC, the value of the complementary
    !    error function.
    !
    implicit none

    integer ( kind = int32 ) :: ind
    real ( kind = dp ) :: x
    real ( kind = dp ) :: error_fc

    real ( kind = dp ) :: ax,bot,e,t,top,w
    real ( kind = dp ), parameter :: c = 0.564189583547756D+00

    real ( kind = dp ), dimension ( 5 ) :: a = (/ &
        0.771058495001320D-04,  -0.133733772997339D-02, &
        0.323076579225834D-01,   0.479137145607681D-01, &
        0.128379167095513D+00 /)
    real ( kind = dp ), dimension(3) :: b = (/ &
        0.301048631703895D-02, &
        0.538971687740286D-01, &
        0.375795757275549D+00 /)
    real ( kind = dp ), dimension ( 8 ) :: p = (/ &
        -1.36864857382717D-07, 5.64195517478974D-01, &
        7.21175825088309D+00, 4.31622272220567D+01, &
        1.52989285046940D+02, 3.39320816734344D+02, &
        4.51918953711873D+02, 3.00459261020162D+02 /)
    real ( kind = dp ), dimension ( 8 ) :: q = (/  &
        1.00000000000000D+00, 1.27827273196294D+01, &
        7.70001529352295D+01, 2.77585444743988D+02, &
        6.38980264465631D+02, 9.31354094850610D+02, &
        7.90950925327898D+02, 3.00459260956983D+02 /)
    real ( kind = dp ), dimension ( 5 ) :: r = (/ &
        2.10144126479064D+00, 2.62370141675169D+01, &
        2.13688200555087D+01, 4.65807828718470D+00, &
        2.82094791773523D-01 /)
    real ( kind = dp ), dimension ( 4 ) :: s = (/ &
        9.41537750555460D+01, 1.87114811799590D+02, &
        9.90191814623914D+01, 1.80124575948747D+02 /)

    !
    !  ABS ( X ) <= 0.5
    !
    ax = abs ( x )

    if ( ax <= 0.5D+00 ) then

        t = x * x

        top = (((( a(1) * t + a(2) ) * t + a(3) ) * t + a(4) ) * t + a(5) ) &
            + 1.0D+00

        bot = (( b(1) * t + b(2) ) * t + b(3) ) * t + 1.0D+00

        error_fc = 0.5D+00 + ( 0.5D+00 &
            - x * ( top / bot ) )

        if ( ind /= 0 ) then
            error_fc = exp ( t ) * error_fc
        end if

        return

    end if
    !
    !  0.5 < abs ( X ) <= 4
    !
    if ( ax <= 4.0D+00 ) then

        top = (((((( p(1) * ax + p(2)) * ax + p(3)) * ax + p(4)) * ax &
            + p(5)) * ax + p(6)) * ax + p(7)) * ax + p(8)

        bot = (((((( q(1) * ax + q(2)) * ax + q(3)) * ax + q(4)) * ax &
            + q(5)) * ax + q(6)) * ax + q(7)) * ax + q(8)

        error_fc = top / bot
        !
        !  4 < ABS ( X )
        !
    else

        if ( x <= -5.6D+00 ) then

            if ( ind == 0 ) then
                error_fc =  2.0D+00
            else
                error_fc =  2.0D+00  * exp ( x * x )
            end if

            return

        end if

        if ( ind == 0 ) then

            if ( 100.0D+00 < x ) then
                error_fc = 0.0D+00
                return
            end if

            if ( -exparg ( 1 ) < x * x ) then
                error_fc = 0.0D+00
                return
            end if

        end if

        t = ( 1.0D+00 / x )**2

        top = ((( r(1) * t + r(2) ) * t + r(3) ) * t + r(4) ) * t + r(5)

        bot = ((( s(1) * t + s(2) ) * t + s(3) ) * t + s(4) ) * t &
            + 1.0D+00

        error_fc = ( c - t * top / bot ) / ax

    end if
    !
    !  Final assembly.
    !
    if ( ind /= 0 ) then

        if ( x < 0.0D+00 ) then
            error_fc =  2.0D+00  * exp ( x * x ) - error_fc
        end if

    else

        w = x * x
        t = w
        e = w - t
        error_fc = (( 0.5D+00 &
            + ( 0.5D+00 - e ) ) * exp ( - t ) ) * error_fc

        if ( x < 0.0D+00 ) then
            error_fc =  2.0D+00  - error_fc
        end if

    end if

    return
    end function error_fc

    function esum ( mu, x )

    !*****************************************************************************80
    !
    !! ESUM evaluates exp ( MU + X ).
    !
    !  Author:
    !
    !    Armido DiDinato, Alfred Morris
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = int32 ) MU, part of the argument.
    !
    !    Input, real ( kind = dp ) X, part of the argument.
    !
    !    Output, real ( kind = dp ) ESUM, the value of exp ( MU + X ).
    !
    implicit none

    integer ( kind = int32 ) :: mu
    real ( kind = dp ) :: x
    real ( kind = dp ) :: esum

    real ( kind = dp ) :: w

    if ( x <= 0.0D+00 ) then
        if ( 0 <= mu ) then
            w = mu + x
            if ( w <= 0.0D+00 ) then
                esum = exp ( w )
                return
            end if
        end if
    else if ( 0.0D+00 < x ) then
        if ( mu <= 0 ) then
            w = mu + x
            if ( 0.0D+00 <= w ) then
                esum = exp ( w )
                return
            end if
        end if
    end if

    w = mu
    esum = exp ( w ) * exp ( x )

    return
    end function esum

    function eval_pol ( a, n, x )

    !*****************************************************************************80
    !
    !! EVAL_POL evaluates a polynomial at X.
    !
    !  Discussion:
    !
    !    EVAL_POL = A(0) + A(1)*X + ... + A(N)*X**N
    !
    !  Modified:
    !
    !    15 December 1999
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) A(0:N), coefficients of the polynomial.
    !
    !    Input, integer ( kind = int32 ) N, length of A.
    !
    !    Input, real ( kind = dp ) X, the point at which the polynomial
    !    is to be evaluated.
    !
    !    Output, real ( kind = dp ) EVAL_POL, the value of the polynomial at X.
    !
    implicit none

    integer ( kind = int32 ) :: n
    real ( kind = dp ) :: a(0:n)
    real ( kind = dp ) :: x
    real ( kind = dp ) :: eval_pol

    real ( kind = dp ) :: term
    integer ( kind = int32 ) :: i

    term = a(n)
    do i = n - 1, 0, -1
        term = term * x + a(i)
    end do

    eval_pol = term

    return
    end function eval_pol

    function exparg ( l )

    !*****************************************************************************80
    !
    !! EXPARG returns the largest or smallest legal argument for EXP.
    !
    !  Discussion:
    !
    !    Only an approximate limit for the argument of EXP is desired.
    !
    !  Modified:
    !
    !    09 December 1999
    !
    !  Author:
    !
    !    Armido DiDinato, Alfred Morris
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = int32 ) L, indicates which limit is desired.
    !    If L = 0, then the largest positive argument for EXP is desired.
    !    Otherwise, the largest negative argument for EXP for which the
    !    result is nonzero is desired.
    !
    !    Output, real ( kind = dp ) EXPARG, the desired value.
    !
    implicit none

    integer ( kind = int32 ) :: l
    real ( kind = dp ) :: exparg

    real ( kind = dp ) :: lnb
    integer ( kind = int32 ) :: b, m

    !
    !  Get the arithmetic base.
    !
    b = ipmpar(4)
    !
    !  Compute the logarithm of the arithmetic base.
    !
    if ( b == 2 ) then
        lnb = 0.69314718055995D+00
    else if ( b == 8 ) then
        lnb = 2.0794415416798D+00
    else if ( b == 16 ) then
        lnb = 2.7725887222398D+00
    else
        lnb = log ( real ( b, kind = dp ) )
    end if

    if ( l /= 0 ) then
        m = ipmpar(9) - 1
        exparg = 0.99999D+00 * ( m * lnb )
    else
        m = ipmpar(10)
        exparg = 0.99999D+00 * ( m * lnb )
    end if

    return
    end function exparg

    subroutine f_cdf_values ( n_data, a, b, x, fx )

    !*****************************************************************************80
    !
    !! F_CDF_VALUES returns some values of the F CDF test function.
    !
    !  Discussion:
    !
    !    The value of F_CDF ( DFN, DFD, X ) can be evaluated in Mathematica by
    !    commands like:
    !
    !      Needs["Statistics`ContinuousDistributions`"]
    !      CDF[FRatioDistribution[ DFN, DFD ], X ]
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    11 June 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Milton Abramowitz, Irene Stegun,
    !    Handbook of Mathematical Functions,
    !    US Department of Commerce, 1964.
    !
    !    Stephen Wolfram,
    !    The Mathematica Book,
    !    Fourth Edition,
    !    Wolfram Media / Cambridge University Press, 1999.
    !
    !  Parameters:
    !
    !    Input/output, integer ( kind = int32 ) N_DATA.  The user sets N_DATA to 0
    !    before the first call.  On each call, the routine increments N_DATA by 1,
    !    and returns the corresponding data; when there is no more data, the
    !    output value of N_DATA will be 0 again.
    !
    !    Output, integer ( kind = int32 ) A, integer B, real ( kind = dp ) X, the
    !    arguments of the function.
    !
    !    Output, real ( kind = dp ) FX, the value of the function.
    !
    implicit none

    integer ( kind = int32 ) :: n_data,a,b
    real ( kind = dp ) :: x,fx

    integer ( kind = int32 ), parameter :: n_max = 20

    integer ( kind = int32 ), save, dimension ( n_max ) :: a_vec = (/ &
        1, 1, 5, 1, &
        2, 4, 1, 6, &
        8, 1, 3, 6, &
        1, 1, 1, 1, &
        2, 3, 4, 5 /)
    integer ( kind = int32 ), save, dimension ( n_max ) :: b_vec = (/ &
        1,  5,  1,  5, &
        10, 20,  5,  6, &
        16,  5, 10, 12, &
        5,  5,  5,  5, &
        5,  5,  5,  5 /)
    real ( kind = dp ), save, dimension ( n_max ) :: fx_vec = (/ &
        0.500000D+00, 0.499971D+00, 0.499603D+00, 0.749699D+00, &
        0.750466D+00, 0.751416D+00, 0.899987D+00, 0.899713D+00, &
        0.900285D+00, 0.950025D+00, 0.950057D+00, 0.950193D+00, &
        0.975013D+00, 0.990002D+00, 0.994998D+00, 0.999000D+00, &
        0.568799D+00, 0.535145D+00, 0.514343D+00, 0.500000D+00 /)
    real ( kind = dp ), save, dimension ( n_max ) :: x_vec = (/ &
        1.00D+00,  0.528D+00, 1.89D+00,  1.69D+00, &
        1.60D+00,  1.47D+00,  4.06D+00,  3.05D+00, &
        2.09D+00,  6.61D+00,  3.71D+00,  3.00D+00, &
        10.01D+00, 16.26D+00, 22.78D+00, 47.18D+00, &
        1.00D+00,  1.00D+00,  1.00D+00,  1.00D+00 /)

    if ( n_data < 0 ) then
        n_data = 0
    end if

    n_data = n_data + 1

    if ( n_max < n_data ) then
        n_data = 0
        a = 0
        b = 0
        x = 0.0D+00
        fx = 0.0D+00
    else
        a = a_vec(n_data)
        b = b_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
    end if

    return
    end subroutine f_cdf_values

    subroutine f_noncentral_cdf_values ( n_data, a, b, lambda, x, fx )

    !*****************************************************************************80
    !
    !! F_NONCENTRAL_CDF_VALUES returns some values of the F CDF test function.
    !
    !  Discussion:
    !
    !    The value of NONCENTRAL_F_CDF ( DFN, DFD, LAMDA, X ) can be evaluated
    !    in Mathematica by commands like:
    !
    !      Needs["Statistics`ContinuousDistributions`"]
    !      CDF[NoncentralFRatioDistribution[ DFN, DFD, LAMBDA ], X ]
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    12 June 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Milton Abramowitz, Irene Stegun,
    !    Handbook of Mathematical Functions,
    !    US Department of Commerce, 1964.
    !
    !    Stephen Wolfram,
    !    The Mathematica Book,
    !    Fourth Edition,
    !    Wolfram Media / Cambridge University Press, 1999.
    !
    !  Parameters:
    !
    !    Input/output, integer ( kind = int32 ) N_DATA.  The user sets N_DATA to 0
    !    before the first call.  On each call, the routine increments N_DATA by 1,
    !    and returns the corresponding data; when there is no more data, the
    !    output value of N_DATA will be 0 again.
    !
    !    Output, integer ( kind = int32 ) A, B, real ( kind = dp ) LAMBDA, the
    !    parameters of the function.
    !
    !    Output, real ( kind = dp ) X, the argument of the function.
    !
    !    Output, real ( kind = dp ) FX, the value of the function.
    !
    implicit none

    integer ( kind = int32 ) :: n_data,a,b
    real ( kind = dp ) :: lambda,x,fx

    integer ( kind = int32 ), parameter :: n_max = 22
    integer ( kind = int32 ), save, dimension ( n_max ) :: a_vec = (/ &
        1,  1,  1,  1, &
        1,  1,  1,  1, &
        1,  1,  2,  2, &
        3,  3,  4,  4, &
        5,  5,  6,  6, &
        8, 16 /)
    integer ( kind = int32 ), save, dimension ( n_max ) :: b_vec = (/ &
        1,  5,  5,  5, &
        5,  5,  5,  5, &
        5,  5,  5, 10, &
        5,  5,  5,  5, &
        1,  5,  6, 12, &
        16,  8 /)
    real ( kind = dp ), save, dimension ( n_max ) :: fx_vec = (/ &
        0.500000D+00, 0.636783D+00, 0.584092D+00, 0.323443D+00, &
        0.450119D+00, 0.607888D+00, 0.705928D+00, 0.772178D+00, &
        0.819105D+00, 0.317035D+00, 0.432722D+00, 0.450270D+00, &
        0.426188D+00, 0.337744D+00, 0.422911D+00, 0.692767D+00, &
        0.363217D+00, 0.421005D+00, 0.426667D+00, 0.446402D+00, &
        0.844589D+00, 0.816368D+00 /)
    real ( kind = dp ), save, dimension ( n_max ) :: lambda_vec = (/ &
        0.00D+00,  0.000D+00, 0.25D+00,  1.00D+00, &
        1.00D+00,  1.00D+00,  1.00D+00,  1.00D+00, &
        1.00D+00,  2.00D+00,  1.00D+00,  1.00D+00, &
        1.00D+00,  2.00D+00,  1.00D+00,  1.00D+00, &
        0.00D+00,  1.00D+00,  1.00D+00,  1.00D+00, &
        1.00D+00,  1.00D+00 /)
    real ( kind = dp ), save, dimension ( n_max ) :: x_vec = (/ &
        1.00D+00,  1.00D+00, 1.00D+00,  0.50D+00, &
        1.00D+00,  2.00D+00, 3.00D+00,  4.00D+00, &
        5.00D+00,  1.00D+00, 1.00D+00,  1.00D+00, &
        1.00D+00,  1.00D+00, 1.00D+00,  2.00D+00, &
        1.00D+00,  1.00D+00, 1.00D+00,  1.00D+00, &
        2.00D+00,  2.00D+00 /)

    if ( n_data < 0 ) then
        n_data = 0
    end if

    n_data = n_data + 1

    if ( n_max < n_data ) then
        n_data = 0
        a = 0
        b = 0
        lambda = 0.0D+00
        x = 0.0D+00
        fx = 0.0D+00
    else
        a = a_vec(n_data)
        b = b_vec(n_data)
        lambda = lambda_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
    end if

    return
    end subroutine f_noncentral_cdf_values

    function fpser ( a, b, x, eps )

    !*****************************************************************************80
    !
    !! FPSER evaluates IX(A,B)(X) for very small B.
    !
    !  Discussion:
    !
    !    This routine is appropriate for use when
    !
    !      B < min ( EPS, EPS * A )
    !
    !    and
    !
    !      X <= 0.5.
    !
    !  Author:
    !
    !    Armido DiDinato, Alfred Morris
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) A, B, parameters of the function.
    !
    !    Input, real ( kind = dp ) X, the point at which the function is to
    !    be evaluated.
    !
    !    Input, real ( kind = dp ) EPS, a tolerance.
    !
    !    Output, real ( kind = dp ) FPSER, the value of IX(A,B)(X).
    !
    implicit none

    real ( kind = dp ) :: a,b,x,eps
    real ( kind = dp ) :: fpser

    real ( kind = dp ) :: an,c,s,t,tol

    fpser = 1.0D+00

    if ( 1.0D-03 * eps < a ) then
        fpser = 0.0D+00
        t = a * log ( x )
        if ( t < exparg ( 1 ) ) then
            return
        end if
        fpser = exp ( t )
    end if
    !
    !  1/B(A,B) = B
    !
    fpser = ( b / a ) * fpser
    tol = eps / a
    an = a + 1.0D+00
    t = x
    s = t / an

    do

        an = an + 1.0D+00
        t = x * t
        c = t / an
        s = s + c

        if ( abs ( c ) <= tol ) then
            exit
        end if

    end do

    fpser = fpser * ( 1.0D+00 + a * s )

    return
    end function fpser

    function gam1 ( a )

    !*****************************************************************************80
    !
    !! GAM1 computes 1 / GAMMA(A+1) - 1 for -0.5 <= A <= 1.5
    !
    !  Author:
    !
    !    Armido DiDinato, Alfred Morris
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) A, forms the argument of the Gamma function.
    !
    !    Output, real ( kind = dp ) GAM1, the value of 1 / GAMMA ( A + 1 ) - 1.
    !
    implicit none

    real ( kind = dp ) :: a
    real ( kind = dp ) :: gam1

    real ( kind = dp ) :: bot,d,t,top,w

    real ( kind = dp ), parameter :: s1 = 0.273076135303957D+00
    real ( kind = dp ), parameter :: s2 = 0.559398236957378D-01

    real ( kind = dp ), parameter, dimension ( 7 ) :: p = (/ &
        0.577215664901533D+00, -0.409078193005776D+00, &
        -0.230975380857675D+00,  0.597275330452234D-01, &
        0.766968181649490D-02, -0.514889771323592D-02, &
        0.589597428611429D-03 /)
    real ( kind = dp ), dimension ( 5 ) :: q = (/ &
        0.100000000000000D+01, 0.427569613095214D+00, &
        0.158451672430138D+00, 0.261132021441447D-01, &
        0.423244297896961D-02 /)
    real ( kind = dp ), dimension ( 9 ) :: r = (/ &
        -0.422784335098468D+00, -0.771330383816272D+00, &
        -0.244757765222226D+00,  0.118378989872749D+00, &
        0.930357293360349D-03, -0.118290993445146D-01, &
        0.223047661158249D-02,  0.266505979058923D-03, &
        -0.132674909766242D-03 /)


    d = a - 0.5D+00

    if ( 0.0D+00 < d ) then
        t = d - 0.5D+00
    else
        t = a
    end if

    if ( t == 0.0D+00 ) then

        gam1 = 0.0D+00

    else if ( 0.0D+00 < t ) then

        top = (((((    &
            p(7)   &
            * t + p(6) ) &
            * t + p(5) ) &
            * t + p(4) ) &
            * t + p(3) ) &
            * t + p(2) ) &
            * t + p(1)

        bot = ((( q(5) * t + q(4) ) * t + q(3) ) * t + q(2) ) * t &
            + 1.0D+00

        w = top / bot

        if ( d <= 0.0D+00 ) then
            gam1 = a * w
        else
            gam1 = ( t / a ) * ( ( w - 0.5D+00 ) &
                - 0.5D+00 )
        end if

    else if ( t < 0.0D+00 ) then

        top = (((((((  &
            r(9)   &
            * t + r(8) ) &
            * t + r(7) ) &
            * t + r(6) ) &
            * t + r(5) ) &
            * t + r(4) ) &
            * t + r(3) ) &
            * t + r(2) ) &
            * t + r(1)

        bot = ( s2 * t + s1 ) * t + 1.0D+00
        w = top / bot

        if ( d <= 0.0D+00 ) then
            gam1 = a * ( ( w + 0.5D+00 ) + 0.5D+00 )
        else
            gam1 = t * w / a
        end if

    end if

    return
    end function gam1

    function gamma_cdflib ( a )

    !*****************************************************************************80
    !
    !! GAMMA_CDFLIB evaluates the gamma function.
    !
    !  Author:
    !
    !    Alfred Morris,
    !    Naval Surface Weapons Center,
    !    Dahlgren, Virginia.
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) A, the argument of the Gamma function.
    !
    !    Output, real ( kind = dp ) GAMMA, the value of the Gamma function.
    !
    implicit none

    real ( kind = dp ) :: a
    real ( kind = dp ) :: gamma_cdflib

    real ( kind = dp ) :: bot,g,lnx,s,t,top,w,x,z
    integer ( kind = int32 ) :: i,j,m,n

    real ( kind = dp ), parameter :: d = 0.41893853320467274178D+00
    real ( kind = dp ), parameter :: pi = 3.1415926535898D+00
    real ( kind = dp ), parameter :: r1 =  0.820756370353826D-03
    real ( kind = dp ), parameter :: r2 = -0.595156336428591D-03
    real ( kind = dp ), parameter :: r3 =  0.793650663183693D-03
    real ( kind = dp ), parameter :: r4 = -0.277777777770481D-02
    real ( kind = dp ), parameter :: r5 =  0.833333333333333D-01

    real ( kind = dp ), dimension ( 7 ) :: p = (/ &
        0.539637273585445D-03, 0.261939260042690D-02, &
        0.204493667594920D-01, 0.730981088720487D-01, &
        0.279648642639792D+00, 0.553413866010467D+00, &
        1.0D+00 /)
    real ( kind = dp ), dimension ( 7 ) :: q = (/ &
        -0.832979206704073D-03,  0.470059485860584D-02, &
        0.225211131035340D-01, -0.170458969313360D+00, &
        -0.567902761974940D-01,  0.113062953091122D+01, &
        1.0D+00 /)

    gamma_cdflib = 0.0D+00
    x = a

    if ( abs ( a ) < 15.0D+00 ) then
        !
        !  Evaluation of GAMMA(A) for |A| < 15
        !
        t = 1.0D+00
        m = int ( a ) - 1
        !
        !  Let T be the product of A-J when 2 <= A.
        !
        if ( 0 <= m ) then

            do j = 1, m
                x = x - 1.0D+00
                t = x * t
            end do

            x = x - 1.0D+00
            !
            !  Let T be the product of A+J WHEN A < 1
            !
        else

            t = a

            if ( a <= 0.0D+00 ) then

                m = - m - 1

                do j = 1, m
                    x = x + 1.0D+00
                    t = x * t
                end do

                x = ( x + 0.5D+00 ) + 0.5D+00
                t = x * t
                if ( t == 0.0D+00 ) then
                    return
                end if

            end if
            !
            !  Check if 1/T can overflow.
            !
            if ( abs ( t ) < 1.0D-30 ) then
                if ( 1.0001D+00 < abs ( t ) * huge ( t ) ) then
                    gamma_cdflib = 1.0D+00 / t
                end if
                return
            end if

        end if
        !
        !  Compute Gamma(1 + X) for 0 <= X < 1.
        !
        top = p(1)
        bot = q(1)
        do i = 2, 7
            top = top * x + p(i)
            bot = bot * x + q(i)
        end do

        gamma_cdflib = top / bot
        !
        !  Termination.
        !
        if ( 1.0D+00 <= a ) then
            gamma_cdflib = gamma_cdflib * t
        else
            gamma_cdflib = gamma_cdflib / t
        end if
        !
        !  Evaluation of Gamma(A) FOR 15 <= ABS ( A ).
        !
    else

        if ( 1000.0D+00 <= abs ( a ) ) then
            return
        end if

        if ( a <= 0.0D+00 ) then

            x = -a
            n = x
            t = x - n

            if ( 0.9D+00 < t ) then
                t = 1.0D+00 - t
            end if

            s = sin ( pi * t ) / pi

            if ( mod ( n, 2 ) == 0 ) then
                s = -s
            end if

            if ( s == 0.0D+00 ) then
                return
            end if

        end if
        !
        !  Compute the modified asymptotic sum.
        !
        t = 1.0D+00 / ( x * x )

        g = (((( r1 * t + r2 ) * t + r3 ) * t + r4 ) * t + r5 ) / x

        lnx = log ( x )
        !
        !  Final assembly.
        !
        z = x
        g = ( d + g ) + ( z - 0.5D+00 ) &
            * ( lnx - 1.0D+00 )
        w = g
        t = g - real ( w, kind = dp )

        if ( 0.99999D+00 * exparg ( 0 ) < w ) then
            return
        end if

        gamma_cdflib = exp ( w )* ( 1.0D+00 + t )

        if ( a < 0.0D+00 ) then
            gamma_cdflib = ( 1.0D+00 / ( gamma_cdflib * s ) ) / x
        end if

    end if

    return
    end function gamma_cdflib

    subroutine gamma_inc ( a, x, ans, qans, ind )

    !*****************************************************************************80
    !
    !! GAMMA_INC evaluates the incomplete gamma ratio functions P(A,X) and Q(A,X).
    !
    !  Modified:
    !
    !    16 April 2005
    !
    !  Author:
    !
    !    Alfred Morris,
    !    Naval Surface Weapons Center,
    !    Dahlgren, Virginia.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) A, X, the arguments of the incomplete
    !    gamma ratio.  A and X must be nonnegative.  A and X cannot
    !    both be zero.
    !
    !    Output, real ( kind = dp ) ANS, QANS.  On normal output,
    !    ANS = P(A,X) and QANS = Q(A,X).  However, ANS is set to 2 if
    !    A or X is negative, or both are 0, or when the answer is
    !    computationally indeterminate because A is extremely large
    !    and X is very close to A.
    !
    !    Input, integer ( kind = int32 ) IND, indicates the accuracy request:
    !    0, as much accuracy as possible.
    !    1, to within 1 unit of the 6-th significant digit,
    !    otherwise, to within 1 unit of the 3rd significant digit.
    !
    !  Local Parameters:
    !
    !     ALOG10 = LN(10)
    !     RT2PIN = 1/SQRT(2*PI)
    !     RTPI   = SQRT(PI)
    !
    implicit none

    real ( kind = dp ) :: a,x,ans,qans
    integer ( kind = int32 ) :: ind

    real ( kind = dp ) :: a2n,a2nm1,acc,am0,amn,an,an0,apn,b2n,b2nm1
    real ( kind = dp ) :: c,c0,c1,c2,c3,c4,c5,c6,cma
    real ( kind = dp ) :: d10,d20,d30,d40,d50,d60,d70
    real ( kind = dp ) :: e,e0,g,h,j,l,r,rta,rtx,s,sum1,t,t1
    real ( kind = dp ) :: tol,twoa,u,w,y,z,x0

    integer ( kind = int32 ) :: i,iop,m,n,n_max

    real ( kind = dp ), parameter :: rt2pin = 0.398942280401433D+00
    real ( kind = dp ), parameter :: rtpi = 1.77245385090552D+00
    real ( kind = dp ), parameter :: alog10 = 2.30258509299405D+00

    real ( kind = dp ), dimension ( 3 ) :: acc0 = (/ &
        5.0D-15, 5.0D-07, 5.0D-04 /)

    real ( kind = dp ) :: big(3),d0(13),d1(12),d2(10),d3(8),d4(6)
    real ( kind = dp ) :: d5(4),d6(2),e00(3),wk(20),x00(3)

    data big(1)/20.0D+00/,big(2)/14.0D+00/,big(3)/10.0D+00/
    data e00(1)/0.25D-03/,e00(2)/0.25D-01/,e00(3)/0.14D+00/
    data x00(1)/31.0D+00/,x00(2)/17.0D+00/,x00(3)/9.7D+00/
    data d0(1)/0.833333333333333D-01/
    data d0(2)/-0.148148148148148D-01/
    data d0(3)/0.115740740740741D-02/,d0(4)/0.352733686067019D-03/
    data d0(5)/-0.178755144032922D-03/,d0(6)/0.391926317852244D-04/
    data d0(7)/-0.218544851067999D-05/,d0(8)/-0.185406221071516D-05/
    data d0(9)/0.829671134095309D-06/,d0(10)/-0.176659527368261D-06/
    data d0(11)/0.670785354340150D-08/,d0(12)/0.102618097842403D-07/
    data d0(13)/-0.438203601845335D-08/
    data d10/-0.185185185185185D-02/,d1(1)/-0.347222222222222D-02/
    data d1(2)/0.264550264550265D-02/,d1(3)/-0.990226337448560D-03/
    data d1(4)/0.205761316872428D-03/,d1(5)/-0.401877572016461D-06/
    data d1(6)/-0.180985503344900D-04/,d1(7)/0.764916091608111D-05/
    data d1(8)/-0.161209008945634D-05/,d1(9)/0.464712780280743D-08/
    data d1(10)/0.137863344691572D-06/,d1(11)/-0.575254560351770D-07/
    data d1(12)/0.119516285997781D-07/
    data d20/0.413359788359788D-02/,d2(1)/-0.268132716049383D-02/
    data d2(2)/0.771604938271605D-03/,d2(3)/0.200938786008230D-05/
    data d2(4)/-0.107366532263652D-03/,d2(5)/0.529234488291201D-04/
    data d2(6)/-0.127606351886187D-04/,d2(7)/0.342357873409614D-07/
    data d2(8)/0.137219573090629D-05/,d2(9)/-0.629899213838006D-06/
    data d2(10)/0.142806142060642D-06/
    data d30/0.649434156378601D-03/,d3(1)/0.229472093621399D-03/
    data d3(2)/-0.469189494395256D-03/,d3(3)/0.267720632062839D-03/
    data d3(4)/-0.756180167188398D-04/,d3(5)/-0.239650511386730D-06/
    data d3(6)/0.110826541153473D-04/,d3(7)/-0.567495282699160D-05/
    data d3(8)/0.142309007324359D-05/
    data d40/-0.861888290916712D-03/,d4(1)/0.784039221720067D-03/
    data d4(2)/-0.299072480303190D-03/,d4(3)/-0.146384525788434D-05/
    data d4(4)/0.664149821546512D-04/,d4(5)/-0.396836504717943D-04/
    data d4(6)/0.113757269706784D-04/
    data d50/-0.336798553366358D-03/,d5(1)/-0.697281375836586D-04/
    data d5(2)/0.277275324495939D-03/,d5(3)/-0.199325705161888D-03/
    data d5(4)/0.679778047793721D-04/
    data d60/0.531307936463992D-03/,d6(1)/-0.592166437353694D-03/
    data d6(2)/0.270878209671804D-03/
    data d70 / 0.344367606892378D-03/

    e = epsilon ( 1.0D+00 )

    if ( a < 0.0D+00 .or. x < 0.0D+00 ) then
        ans = 2.0D+00
        return
    end if

    if ( a == 0.0D+00 .and. x == 0.0D+00 ) then
        ans = 2.0D+00
        return
    end if

    if ( a * x == 0.0D+00 ) then
        if ( x <= a ) then
            ans = 0.0D+00
            qans = 1.0D+00
        else
            ans = 1.0D+00
            qans = 0.0D+00
        end if
        return
    end if

    iop = ind + 1
    if ( iop /= 1 .and. iop /= 2 ) iop = 3
    acc = max ( acc0(iop), e )
    e0 = e00(iop)
    x0 = x00(iop)
    !
    !  Select the appropriate algorithm.
    !
    if ( 1.0D+00 <= a ) then
        go to 10
    end if

    if ( a == 0.5D+00 ) then
        go to 390
    end if

    if ( x < 1.1D+00 ) then
        go to 160
    end if

    t1 = a * log ( x ) - x
    u = a * exp ( t1 )

    if ( u == 0.0D+00 ) then
        ans = 1.0D+00
        qans = 0.0D+00
        return
    end if

    r = u * ( 1.0D+00 + gam1 ( a ) )
    go to 250

10  continue

    if ( big(iop) <= a ) then
        go to 30
    end if

    if ( x < a .or. x0 <= x ) then
        go to 20
    end if

    twoa = a + a
    m = int ( twoa )

    if ( twoa == real ( m, kind = dp ) ) then
        i = m / 2
        if ( a == real ( i, kind = dp ) ) then
            go to 210
        end if
        go to 220
    end if

20  continue

    t1 = a * log ( x ) - x
    r = exp ( t1 ) / gamma_cdflib ( a )
    go to 40

30  continue

    l = x / a

    if ( l == 0.0D+00 ) then
        ans = 0.0D+00
        qans = 1.0D+00
        return
    end if

    s = 0.5D+00 + ( 0.5D+00 - l )
    z = rlog ( l )
    if ( 700.0D+00 / a <= z ) then
        go to 410
    end if

    y = a * z
    rta = sqrt ( a )

    if ( abs ( s ) <= e0 / rta ) then
        go to 330
    end if

    if ( abs ( s ) <= 0.4D+00 ) then
        go to 270
    end if

    t = ( 1.0D+00 / a )**2
    t1 = ((( 0.75D+00 * t - 1.0D+00 ) * t + 3.5D+00 ) &
        * t - 105.0D+00 ) / ( a * 1260.0D+00 )
    t1 = t1 - y
    r = rt2pin * rta * exp ( t1 )

40  continue

    if ( r == 0.0D+00 ) then
        if ( x <= a ) then
            ans = 0.0D+00
            qans = 1.0D+00
        else
            ans = 1.0D+00
            qans = 0.0D+00
        end if
        return
    end if

    if ( x <= max ( a, alog10 ) ) then
        go to 50
    end if

    if ( x < x0 ) then
        go to 250
    end if

    go to 100
    !
    !  Taylor series for P/R.
    !
50  continue

    apn = a + 1.0D+00
    t = x / apn
    wk(1) = t

    n = 20

    do i = 2, 20
        apn = apn + 1.0D+00
        t = t * ( x / apn )
        if ( t <= 1.0D-03 ) then
            n = i
            exit
        end if
        wk(i) = t
    end do

    sum1 = t

    tol = 0.5D+00 * acc

    do

        apn = apn + 1.0D+00
        t = t * ( x / apn )
        sum1 = sum1 + t

        if ( t <= tol ) then
            exit
        end if

    end do

    n_max = n - 1
    do m = 1, n_max
        n = n - 1
        sum1 = sum1 + wk(n)
    end do

    ans = ( r / a ) * ( 1.0D+00 + sum1 )
    qans = 0.5D+00 + ( 0.5D+00 - ans )
    return
    !
    !  Asymptotic expansion.
    !
100 continue

    amn = a - 1.0D+00
    t = amn / x
    wk(1) = t

    n = 20

    do i = 2, 20
        amn = amn - 1.0D+00
        t = t * ( amn / x )
        if ( abs ( t ) <= 1.0D-03 ) then
            n = i
            exit
        end if
        wk(i) = t
    end do

    sum1 = t

    do

        if ( abs ( t ) <= acc ) then
            exit
        end if

        amn = amn - 1.0D+00
        t = t * ( amn / x )
        sum1 = sum1 + t

    end do

    n_max = n - 1
    do m = 1, n_max
        n = n - 1
        sum1 = sum1 + wk(n)
    end do
    qans = ( r / x ) * ( 1.0D+00 + sum1 )
    ans = 0.5D+00 + ( 0.5D+00 - qans )
    return
    !
    !  Taylor series for P(A,X)/X**A
    !
160 continue

    an = 3.0D+00
    c = x
    sum1 = x / ( a + 3.0D+00 )
    tol = 3.0D+00 * acc / ( a + 1.0D+00 )

    do

        an = an + 1.0D+00
        c = -c * ( x / an )
        t = c / ( a + an )
        sum1 = sum1 + t

        if ( abs ( t ) <= tol ) then
            exit
        end if

    end do

    j = a * x * ( ( sum1 / 6.0D+00 - 0.5D+00 / &
        ( a +  2.0D+00  ) ) * x + 1.0D+00 &
        / ( a + 1.0D+00 ) )

    z = a * log ( x )
    h = gam1 ( a )
    g = 1.0D+00 + h

    if ( x < 0.25D+00 ) then
        go to 180
    end if

    if ( a < x / 2.59D+00 ) then
        go to 200
    end if

    go to 190

180 continue

    if ( -0.13394D+00 < z ) then
        go to 200
    end if

190 continue

    w = exp ( z )
    ans = w * g * ( 0.5D+00 + ( 0.5D+00 - j ))
    qans = 0.5D+00 + ( 0.5D+00 - ans )
    return

200 continue

    l = rexp ( z )
    w = 0.5D+00 + ( 0.5D+00 + l )
    qans = ( w * j - l ) * g - h

    if ( qans < 0.0D+00 ) then
        ans = 1.0D+00
        qans = 0.0D+00
        return
    end if

    ans = 0.5D+00 + ( 0.5D+00 - qans )
    return
    !
    !  Finite sums for Q when 1 <= A and 2*A is an integer.
    !
210 continue

    sum1 = exp ( - x )
    t = sum1
    n = 1
    c = 0.0D+00
    go to 230

220 continue

    rtx = sqrt ( x )
    sum1 = error_fc ( 0, rtx )
    t = exp ( -x ) / ( rtpi * rtx )
    n = 0
    c = -0.5D+00

230 continue

    do while ( n /= i )
        n = n + 1
        c = c + 1.0D+00
        t = ( x * t ) / c
        sum1 = sum1 + t
    end do

240 continue

    qans = sum1
    ans = 0.5D+00 + ( 0.5D+00 - qans )
    return
    !
    !  Continued fraction expansion.
    !
250 continue

    tol = max ( 5.0D+00 * e, acc )
    a2nm1 = 1.0D+00
    a2n = 1.0D+00
    b2nm1 = x
    b2n = x + ( 1.0D+00 - a )
    c = 1.0D+00

    do

        a2nm1 = x * a2n + c * a2nm1
        b2nm1 = x * b2n + c * b2nm1
        am0 = a2nm1 / b2nm1
        c = c + 1.0D+00
        cma = c - a
        a2n = a2nm1 + cma * a2n
        b2n = b2nm1 + cma * b2n
        an0 = a2n / b2n

        if ( abs ( an0 - am0 ) < tol * an0 ) then
            exit
        end if

    end do

    qans = r * an0
    ans = 0.5D+00 + ( 0.5D+00 - qans )
    return
    !
    !  General Temme expansion.
    !
270 continue

    if ( abs ( s ) <= 2.0D+00 * e .and. 3.28D-03 < a * e * e ) then
        ans =  2.0D+00
        return
    end if

    c = exp ( - y )
    w = 0.5D+00 * error_fc ( 1, sqrt ( y ) )
    u = 1.0D+00 / a
    z = sqrt ( z + z )

    if ( l < 1.0D+00 ) then
        z = -z
    end if

    if ( iop < 2 ) then

        if ( abs ( s ) <= 1.0D-03 ) then

            c0 = ((((((     &
                d0(7)   &
                * z + d0(6) ) &
                * z + d0(5) ) &
                * z + d0(4) ) &
                * z + d0(3) ) &
                * z + d0(2) ) &
                * z + d0(1) ) &
                * z - 1.0D+00 / 3.0D+00

            c1 = (((((      &
                d1(6)   &
                * z + d1(5) ) &
                * z + d1(4) ) &
                * z + d1(3) ) &
                * z + d1(2) ) &
                * z + d1(1) ) &
                * z + d10

            c2 = ((((d2(5)*z+d2(4))*z+d2(3))*z+d2(2))*z+d2(1))*z + d20

            c3 = (((d3(4)*z+d3(3))*z+d3(2))*z+d3(1))*z + d30

            c4 = ( d4(2) * z + d4(1) ) * z + d40
            c5 = ( d5(2) * z + d5(1) ) * z + d50
            c6 = d6(1) * z + d60

            t = (((((( d70 &
                * u + c6 ) &
                * u + c5 ) &
                * u + c4 ) &
                * u + c3 ) &
                * u + c2 ) &
                * u + c1 ) &
                * u + c0

        else

            c0 = (((((((((((( &
                d0(13)   &
                * z + d0(12) ) &
                * z + d0(11) ) &
                * z + d0(10) ) &
                * z + d0(9)  ) &
                * z + d0(8)  ) &
                * z + d0(7)  ) &
                * z + d0(6)  ) &
                * z + d0(5)  ) &
                * z + d0(4)  ) &
                * z + d0(3)  ) &
                * z + d0(2)  ) &
                * z + d0(1)  ) &
                * z - 1.0D+00 / 3.0D+00

            c1 = ((((((((((( &
                d1(12) &
                * z + d1(11) &
                ) * z + d1(10) &
                ) * z + d1(9)  &
                ) * z + d1(8)  &
                ) * z + d1(7)  &
                ) * z + d1(6)  &
                ) * z + d1(5)  &
                ) * z + d1(4)  &
                ) * z + d1(3)  &
                ) * z + d1(2)  &
                ) * z + d1(1)  &
                ) * z + d10

            c2 = ((((((((( &
                d2(10) &
                * z + d2(9) &
                ) * z + d2(8) &
                ) * z + d2(7) &
                ) * z + d2(6) &
                ) * z + d2(5) &
                ) * z + d2(4) &
                ) * z + d2(3) &
                ) * z + d2(2) &
                ) * z + d2(1) &
                ) * z + d20

            c3 = ((((((( &
                d3(8) &
                * z + d3(7) &
                ) * z + d3(6) &
                ) * z + d3(5) &
                ) * z + d3(4) &
                ) * z + d3(3) &
                ) * z + d3(2) &
                ) * z + d3(1) &
                ) * z + d30

            c4 = ((((( d4(6)*z+d4(5))*z+d4(4))*z+d4(3))*z+d4(2))*z+d4(1))*z + d40

            c5 = (((d5(4)*z+d5(3))*z+d5(2))*z+d5(1))*z + d50

            c6 = ( d6(2) * z + d6(1) ) * z + d60

            t = ((((((   &
                d70    &
                * u + c6 ) &
                * u + c5 ) &
                * u + c4 ) &
                * u + c3 ) &
                * u + c2 ) &
                * u + c1 ) &
                * u + c0

        end if

    else if ( iop == 2 ) then

        c0 = (((((      &
            d0(6)   &
            * z + d0(5) ) &
            * z + d0(4) ) &
            * z + d0(3) ) &
            * z + d0(2) ) &
            * z + d0(1) ) &
            * z - 1.0D+00 / 3.0D+00

        c1 = ((( d1(4) * z + d1(3) ) * z + d1(2) ) * z + d1(1) ) * z + d10
        c2 = d2(1) * z + d20
        t = ( c2 * u + c1 ) * u + c0

    else if ( 2 < iop ) then

        t = (( d0(3) * z + d0(2) ) * z + d0(1) ) * z - 1.0D+00 / 3.0D+00

    end if

310 continue

    if ( 1.0D+00 <= l ) then
        qans = c * ( w + rt2pin * t / rta )
        ans = 0.5D+00 + ( 0.5D+00 - qans )
    else
        ans = c * ( w - rt2pin * t / rta )
        qans = 0.5D+00 + ( 0.5D+00 - ans )
    end if

    return
    !
    !  Temme expansion for L = 1
    !
330 continue

    if ( 3.28D-03 < a * e * e ) then
        ans =  2.0D+00
        return
    end if

    c = 0.5D+00 + ( 0.5D+00 - y )
    w = ( 0.5D+00 - sqrt ( y ) &
        * ( 0.5D+00 &
        + ( 0.5D+00 - y / 3.0D+00 ) ) / rtpi ) / c
    u = 1.0D+00 / a
    z = sqrt ( z + z )

    if ( l < 1.0D+00 ) then
        z = -z
    end if

    if ( iop < 2 ) then

        c0 = ((((((     &
            d0(7)   &
            * z + d0(6) ) &
            * z + d0(5) ) &
            * z + d0(4) ) &
            * z + d0(3) ) &
            * z + d0(2) ) &
            * z + d0(1) ) &
            * z - 1.0D+00 / 3.0D+00

        c1 = (((((      &
            d1(6)   &
            * z + d1(5) ) &
            * z + d1(4) ) &
            * z + d1(3) ) &
            * z + d1(2) ) &
            * z + d1(1) ) &
            * z + d10

        c2 = ((((d2(5)*z+d2(4))*z+d2(3))*z+d2(2))*z+d2(1))*z + d20

        c3 = (((d3(4)*z+d3(3))*z+d3(2))*z+d3(1))*z + d30

        c4 = ( d4(2) * z + d4(1) ) * z + d40
        c5 = ( d5(2) * z + d5(1) ) * z + d50
        c6 = d6(1) * z + d60

        t = (((((( d70 &
            * u + c6 ) &
            * u + c5 ) &
            * u + c4 ) &
            * u + c3 ) &
            * u + c2 ) &
            * u + c1 ) &
            * u + c0

    else if ( iop == 2 ) then

        c0 = ( d0(2) * z + d0(1) ) * z - 1.0D+00 / 3.0D+00
        c1 = d1(1) * z + d10
        t = ( d20 * u + c1 ) * u + c0

    else if ( 2 < iop ) then

        t = d0(1) * z - 1.0D+00 / 3.0D+00

    end if

    go to 310
    !
    !  Special cases
    !
390 continue

    if ( x < 0.25D+00 ) then
        ans = error_f ( sqrt ( x ) )
        qans = 0.5D+00 + ( 0.5D+00 - ans )
    else
        qans = error_fc ( 0, sqrt ( x ) )
        ans = 0.5D+00 + ( 0.5D+00 - qans )
    end if

    return

410 continue

    if ( abs ( s ) <= 2.0D+00 * e ) then
        ans =  2.0D+00
        return
    end if

    if ( x <= a ) then
        ans = 0.0D+00
        qans = 1.0D+00
    else
        ans = 1.0D+00
        qans = 0.0D+00
    end if

    return
    end subroutine gamma_inc

    subroutine gamma_inc_inv ( a, x, x0, p, q, ierr )

    !*****************************************************************************80
    !
    !! GAMMA_INC_INV computes the inverse incomplete gamma ratio function.
    !
    !  Discussion:
    !
    !    The routine is given positive A, and nonnegative P and Q where P + Q = 1.
    !    The value X is computed with the property that P(A,X) = P and Q(A,X) = Q.
    !    Schroder iteration is employed.  The routine attempts to compute X
    !    to 10 significant digits if this is possible for the particular computer
    !    arithmetic being used.
    !
    !  Author:
    !
    !    Alfred Morris,
    !    Naval Surface Weapons Center,
    !    Dahlgren, Virginia.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) A, the parameter in the incomplete gamma
    !    ratio.  A must be positive.
    !
    !    Output, real ( kind = dp ) X, the computed point for which the
    !    incomplete gamma functions have the values P and Q.
    !
    !    Input, real ( kind = dp ) X0, an optional initial approximation
    !    for the solution X.  If the user does not want to supply an
    !    initial approximation, then X0 should be set to 0, or a negative
    !    value.
    !
    !    Input, real ( kind = dp ) P, Q, the values of the incomplete gamma
    !    functions, for which the corresponding argument is desired.
    !
    !    Output, integer ( kind = int32 ) IERR, error flag.
    !    0, the solution was obtained. Iteration was not used.
    !    0 < K, The solution was obtained. IERR iterations were performed.
    !    -2, A <= 0
    !    -3, No solution was obtained. The ratio Q/A is too large.
    !    -4, P + Q /= 1
    !    -6, 20 iterations were performed. The most recent value obtained
    !        for X is given.  This cannot occur if X0 <= 0.
    !    -7, Iteration failed. No value is given for X.
    !        This may occur when X is approximately 0.
    !    -8, A value for X has been obtained, but the routine is not certain
    !        of its accuracy.  Iteration cannot be performed in this
    !        case. If X0 <= 0, this can occur only when P or Q is
    !        approximately 0. If X0 is positive then this can occur when A is
    !        exceedingly close to X and A is extremely large (say A .GE. 1.E20).
    !
    implicit none

    real ( kind = dp ) :: a,x,x0,p,q
    integer ( kind = int32 ) :: ierr

    real ( kind = dp ) :: am1,amax,ap1,ap2,ap3,apn,b
    real ( kind = dp ) :: c1,c2,c3,c4,c5,d,e,e2,eps
    real ( kind = dp ) :: g,h,pn,qg,qn,r,rta,s,s2,sum1
    real ( kind = dp ) :: t,u,w,xn,y,z
    integer ( kind = int32 ) :: iop

    real ( kind = dp ), parameter :: a0 = 3.31125922108741D+00
    real ( kind = dp ), parameter :: a1 = 11.6616720288968D+00
    real ( kind = dp ), parameter :: a2 = 4.28342155967104D+00
    real ( kind = dp ), parameter :: a3 = 0.213623493715853D+00
    real ( kind = dp ), parameter :: b1 = 6.61053765625462D+00
    real ( kind = dp ), parameter :: b2 = 6.40691597760039D+00
    real ( kind = dp ), parameter :: b3 = 1.27364489782223D+00
    real ( kind = dp ), parameter :: b4 = .036117081018842D+00
    real ( kind = dp ), parameter :: c = 0.577215664901533D+00
    real ( kind = dp ), parameter :: tol = 1.0D-05
    real ( kind = dp ), parameter :: two =  2.0D+00
    real ( kind = dp ), parameter :: half = 0.5D+00
    real ( kind = dp ), parameter :: ln10 = 2.302585D+00

    real ( kind = dp ), dimension(2) :: amin = (/ &
        500.0D+00, 100.0D+00 /)
    real ( kind = dp ), dimension ( 2 ) :: bmin = (/ &
        1.0D-28, 1.0D-13 /)
    real ( kind = dp ), dimension ( 2 ) :: dmin = (/ &
        1.0D-06, 1.0D-04 /)
    real ( kind = dp ), dimension ( 2 ) :: emin = (/ &
        2.0D-03, 6.0D-03 /)
    real ( kind = dp ), dimension ( 2 ) :: eps0 = (/ &
        1.0D-10, 1.0D-08 /)

    e = epsilon ( e )

    x = 0.0D+00

    if ( a <= 0.0D+00 ) then
        ierr = -2
        return
    end if

    t = p + q - 1.0D+00

    if ( e < abs ( t ) ) then
        ierr = -4
        return
    end if

    ierr = 0

    if ( p == 0.0D+00 ) then
        return
    end if

    if ( q == 0.0D+00 ) then
        x = huge ( x )
        return
    end if

    if ( a == 1.0D+00 ) then
        if ( 0.9D+00 <= q ) then
            x = -alnrel ( - p )
        else
            x = -log ( q )
        end if
        return
    end if

    e2 = two * e
    amax = 0.4D-10 / ( e * e )

    if ( 1.0D-10 < e ) then
        iop = 2
    else
        iop = 1
    end if

    eps = eps0(iop)
    xn = x0

    if ( 0.0D+00 < x0 ) then
        go to 160
    end if
    !
    !  Selection of the initial approximation XN of X when A < 1.
    !
    if ( 1.0D+00 < a ) then
        go to 80
    end if

    g = gamma_cdflib ( a + 1.0D+00 )
    qg = q * g

    if ( qg == 0.0D+00 ) then
        x = huge ( x )
        ierr = -8
        return
    end if

    b = qg / a

    if ( 0.6D+00 * a < qg ) then
        go to 40
    end if

    if ( a < 0.30D+00 .and. 0.35D+00 <= b ) then
        t = exp ( - ( b + c ) )
        u = t * exp ( t )
        xn = t * exp ( u )
        go to 160
    end if

    if ( 0.45D+00 <= b ) then
        go to 40
    end if

    if ( b == 0.0D+00 ) then
        x = huge ( x )
        ierr = -8
        return
    end if

    y = -log ( b )
    s = half + ( half - a )
    z = log ( y )
    t = y - s * z

    if ( 0.15D+00 <= b ) then
        xn = y - s * log ( t ) - log ( 1.0D+00 + s / ( t + 1.0D+00 ) )
        go to 220
    end if

    if ( 0.01D+00 < b ) then
        u = ( ( t + two * ( 3.0D+00 - a ) ) * t &
            + ( two - a ) * ( 3.0D+00 - a )) / &
            ( ( t + ( 5.0D+00 - a ) ) * t + two )
        xn = y - s * log ( t ) - log ( u )
        go to 220
    end if

30  continue

    c1 = -s * z
    c2 = -s * ( 1.0D+00 + c1 )

    c3 = s * (( half * c1 &
        + ( two - a ) ) * c1 + ( 2.5D+00 - 1.5D+00 * a ) )

    c4 = -s * ((( c1 / 3.0D+00 + ( 2.5D+00 - 1.5D+00 * a ) ) * c1 &
        + ( ( a - 6.0D+00 ) * a + 7.0D+00 ) ) &
        * c1 + ( ( 11.0D+00 * a - 46.0D+00 ) * a + 47.0D+00 ) / 6.0D+00 )

    c5 = -s * (((( - c1 / 4.0D+00 + ( 11.0D+00 * a - 17.0D+00 ) / 6.0D+00 ) * c1 &
        + ( ( -3.0D+00 * a + 13.0D+00 ) * a - 13.0D+00 ) ) * c1 &
        + half &
        * ( ( ( two * a - 25.0D+00 ) * a + 72.0D+00 ) &
        * a - 61.0D+00 ) ) * c1 &
        + ( ( ( 25.0D+00 * a - 195.0D+00 ) * a &
        + 477.0D+00 ) * a - 379.0D+00 ) / 12.0D+00 )

    xn = (((( c5 / y + c4 ) / y + c3 ) / y + c2 ) / y + c1 ) + y

    if ( 1.0D+00 < a ) then
        go to 220
    end if

    if ( bmin(iop) < b ) then
        go to 220
    end if

    x = xn
    return

40  continue

    if ( b * q <= 1.0D-08 ) then
        xn = exp ( - ( q / a + c ))
    else if ( 0.9D+00 < p ) then
        xn = exp ( ( alnrel ( - q ) + gamma_ln1 ( a )) / a )
    else
        xn = exp ( log ( p * g ) / a )
    end if

    if ( xn == 0.0D+00 ) then
        ierr = -3
        return
    end if

    t = half + ( half - xn / ( a + 1.0D+00 ))
    xn = xn / t
    go to 160
    !
    !  Selection of the initial approximation XN of X when 1 < A.
    !
80  continue

    if ( 0.5D+00 < q ) then
        w = log ( p )
    else
        w = log ( q )
    end if

    t = sqrt ( - two * w )

    s = t - ((( a3 * t + a2 ) * t + a1 ) * t + a0 ) / (((( &
        b4 * t + b3 ) * t + b2 ) * t + b1 ) * t + 1.0D+00 )

    if ( 0.5D+00 < q ) then
        s = -s
    end if

    rta = sqrt ( a )
    s2 = s * s

    xn = a + s * rta + ( s2 - 1.0D+00 ) / 3.0D+00 + s * ( s2 - 7.0D+00 ) &
        / ( 36.0D+00 * rta ) - ( ( 3.0D+00 * s2 + 7.0D+00 ) * s2 - 16.0D+00 ) &
        / ( 810.0D+00 * a ) + s * (( 9.0D+00 * s2 + 256.0D+00 ) * s2 - 433.0D+00 ) &
        / ( 38880.0D+00 * a * rta )

    xn = max ( xn, 0.0D+00 )

    if ( amin(iop) <= a ) then

        x = xn
        d = half + ( half - x / a )

        if ( abs ( d ) <= dmin(iop) ) then
            return
        end if

    end if

110 continue

    if ( p <= 0.5D+00 ) then
        go to 130
    end if

    if ( xn < 3.0D+00 * a ) then
        go to 220
    end if

    y = - ( w + gamma_log ( a ) )
    d = max ( two, a * ( a - 1.0D+00 ) )

    if ( ln10 * d <= y ) then
        s = 1.0D+00 - a
        z = log ( y )
        go to 30
    end if

120 continue

    t = a - 1.0D+00
    xn = y + t * log ( xn ) - alnrel ( -t / ( xn + 1.0D+00 ) )
    xn = y + t * log ( xn ) - alnrel ( -t / ( xn + 1.0D+00 ) )
    go to 220

130 continue

    ap1 = a + 1.0D+00

    if ( 0.70D+00 * ap1 < xn ) then
        go to 170
    end if

    w = w + gamma_log ( ap1 )

    if ( xn <= 0.15 * ap1 ) then
        ap2 = a + two
        ap3 = a + 3.0D+00
        x = exp ( ( w + x ) / a )
        x = exp ( ( w + x - log ( 1.0D+00 + ( x / ap1 ) &
            * ( 1.0D+00 + x / ap2 ) ) ) / a )
        x = exp ( ( w + x - log ( 1.0D+00 + ( x / ap1 ) &
            * ( 1.0D+00 + x / ap2 ) ) ) / a )
        x = exp ( ( w + x - log ( 1.0D+00 + ( x / ap1 ) &
            * ( 1.0D+00 + ( x / ap2 ) &
            * ( 1.0D+00 + x / ap3 ) ) ) ) / a )
        xn = x

        if ( xn <= 1.0D-02 * ap1 ) then
            if ( xn <= emin(iop) * ap1 ) then
                return
            end if
            go to 170
        end if

    end if

    apn = ap1
    t = xn / apn
    sum1 = 1.0D+00 + t

    do

        apn = apn + 1.0D+00
        t = t * ( xn / apn )
        sum1 = sum1 + t

        if ( t <= 1.0D-04 ) then
            exit
        end if

    end do

    t = w - log ( sum1 )
    xn = exp ( ( xn + t ) / a )
    xn = xn * ( 1.0D+00 - ( a * log ( xn ) - xn - t ) / ( a - xn ) )
    go to 170
    !
    !  Schroder iteration using P.
    !
160 continue

    if ( 0.5D+00 < p ) then
        go to 220
    end if

170 continue

    if ( p <= 1.0D+10 * tiny ( p ) ) then
        x = xn
        ierr = -8
        return
    end if

    am1 = ( a - half ) - half

180 continue

    if ( amax < a ) then
        d = half + ( half - xn / a )
        if ( abs ( d ) <= e2 ) then
            x = xn
            ierr = -8
            return
        end if
    end if

190 continue

    if ( 20 <= ierr ) then
        ierr = -6
        return
    end if

    ierr = ierr + 1
    call gamma_inc ( a, xn, pn, qn, 0 )

    if ( pn == 0.0D+00 .or. qn == 0.0D+00 ) then
        x = xn
        ierr = -8
        return
    end if

    r = rcomp ( a, xn )

    if ( r == 0.0D+00 ) then
        x = xn
        ierr = -8
        return
    end if

    t = ( pn - p ) / r
    w = half * ( am1 - xn )

    if ( abs ( t ) <= 0.1D+00 .and. abs ( w * t ) <= 0.1D+00 ) then
        go to 200
    end if

    x = xn * ( 1.0D+00 - t )

    if ( x <= 0.0D+00 ) then
        ierr = -7
        return
    end if

    d = abs ( t )
    go to 210

200 continue

    h = t * ( 1.0D+00 + w * t )
    x = xn * ( 1.0D+00 - h )

    if ( x <= 0.0D+00 ) then
        ierr = -7
        return
    end if

    if ( 1.0D+00 <= abs ( w ) .and. abs ( w ) * t * t <= eps ) then
        return
    end if

    d = abs ( h )

210 continue

    xn = x

    if ( d <= tol ) then

        if ( d <= eps ) then
            return
        end if

        if ( abs ( p - pn ) <= tol * p ) then
            return
        end if

    end if

    go to 180
    !
    !  Schroder iteration using Q.
    !
220 continue

    if ( q <= 1.0D+10 * tiny ( q ) ) then
        x = xn
        ierr = -8
        return
    end if

    am1 = ( a - half ) - half

230 continue

    if ( amax < a ) then
        d = half + ( half - xn / a )
        if ( abs ( d ) <= e2 ) then
            x = xn
            ierr = -8
            return
        end if
    end if

    if ( 20 <= ierr ) then
        ierr = -6
        return
    end if

    ierr = ierr + 1
    call gamma_inc ( a, xn, pn, qn, 0 )

    if ( pn == 0.0D+00 .or. qn == 0.0D+00 ) then
        x = xn
        ierr = -8
        return
    end if

    r = rcomp ( a, xn )

    if ( r == 0.0D+00 ) then
        x = xn
        ierr = -8
        return
    end if

    t = ( q - qn ) / r
    w = half * ( am1 - xn )

    if ( abs ( t ) <= 0.1 .and. abs ( w * t ) <= 0.1 ) then
        go to 250
    end if

    x = xn * ( 1.0D+00 - t )

    if ( x <= 0.0D+00 ) then
        ierr = -7
        return
    end if

    d = abs ( t )
    go to 260

250 continue

    h = t * ( 1.0D+00 + w * t )
    x = xn * ( 1.0D+00 - h )

    if ( x <= 0.0D+00 ) then
        ierr = -7
        return
    end if

    if ( 1.0D+00 <= abs ( w ) .and. abs ( w ) * t * t <= eps ) then
        return
    end if

    d = abs ( h )

260 continue

    xn = x

    if ( tol < d ) then
        go to 230
    end if

    if ( d <= eps ) then
        return
    end if

    if ( abs ( q - qn ) <= tol * q ) then
        return
    end if

    go to 230
    end subroutine gamma_inc_inv

    subroutine gamma_inc_values ( n_data, a, x, fx )

    !*****************************************************************************80
    !
    !! GAMMA_INC_VALUES returns some values of the incomplete Gamma function.
    !
    !  Discussion:
    !
    !    The (normalized) incomplete Gamma function P(A,X) is defined as:
    !
    !      PN(A,X) = 1/GAMMA(A) * Integral ( 0 <= T <= X ) T**(A-1) * exp(-T) dT.
    !
    !    With this definition, for all A and X,
    !
    !      0 <= PN(A,X) <= 1
    !
    !    and
    !
    !      PN(A,INFINITY) = 1.0
    !
    !    Mathematica can compute this value as
    !
    !      1 - GammaRegularized[A,X]
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    08 May 2001
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Milton Abramowitz, Irene Stegun,
    !    Handbook of Mathematical Functions,
    !    US Department of Commerce, 1964.
    !
    !  Parameters:
    !
    !    Input/output, integer ( kind = int32 ) N_DATA.  The user sets N_DATA to 0
    !    before the first call.  On each call, the routine increments N_DATA by 1,
    !    and returns the corresponding data; when there is no more data, the
    !    output value of N_DATA will be 0 again.
    !
    !    Output, real ( kind = dp ) A, X, the arguments of the function.
    !
    !    Output, real ( kind = dp ) FX, the value of the function.
    !
    implicit none

    integer ( kind = int32 ) :: n_data
    real ( kind = dp ) :: a,x,fx

    integer ( kind = int32 ), parameter :: n_max = 20

    real ( kind = dp ), save, dimension ( n_max ) :: a_vec = (/ &
        0.1D+00,  0.1D+00,  0.1D+00,  0.5D+00, &
        0.5D+00,  0.5D+00,  1.0D+00,  1.0D+00, &
        1.0D+00,  1.1D+00,  1.1D+00,  1.1D+00, &
        2.0D+00,  2.0D+00,  2.0D+00,  6.0D+00, &
        6.0D+00, 11.0D+00, 26.0D+00, 41.0D+00 /)
    real ( kind = dp ), save, dimension ( n_max ) :: fx_vec = (/ &
        0.7420263D+00, 0.9119753D+00, 0.9898955D+00, 0.2931279D+00, &
        0.7656418D+00, 0.9921661D+00, 0.0951626D+00, 0.6321206D+00, &
        0.9932621D+00, 0.0757471D+00, 0.6076457D+00, 0.9933425D+00, &
        0.0091054D+00, 0.4130643D+00, 0.9931450D+00, 0.0387318D+00, &
        0.9825937D+00, 0.9404267D+00, 0.4863866D+00, 0.7359709D+00 /)
    real ( kind = dp ), save, dimension ( n_max ) :: x_vec = (/ &
        3.1622777D-02, 3.1622777D-01, 1.5811388D+00, 7.0710678D-02, &
        7.0710678D-01, 3.5355339D+00, 0.1000000D+00, 1.0000000D+00, &
        5.0000000D+00, 1.0488088D-01, 1.0488088D+00, 5.2440442D+00, &
        1.4142136D-01, 1.4142136D+00, 7.0710678D+00, 2.4494897D+00, &
        1.2247449D+01, 1.6583124D+01, 2.5495098D+01, 4.4821870D+01 /)

    if ( n_data < 0 ) then
        n_data = 0
    end if

    n_data = n_data + 1

    if ( n_max < n_data ) then
        n_data = 0
        a = 0.0D+00
        x = 0.0D+00
        fx = 0.0D+00
    else
        a = a_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
    end if

    return
    end subroutine gamma_inc_values

    function gamma_ln1 ( a )

    !*****************************************************************************80
    !
    !! GAMMA_LN1 evaluates ln ( Gamma ( 1 + A ) ), for -0.2 <= A <= 1.25.
    !
    !  Author:
    !
    !    Armido DiDinato, Alfred Morris
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) A, defines the argument of the function.
    !
    !    Output, real ( kind = dp ) GAMMA_LN1, the value of ln ( Gamma ( 1 + A ) ).
    !
    implicit none

    real ( kind = dp ) :: a
    real ( kind = dp ) :: gamma_ln1

    real ( kind = dp ) :: bot,top,x

    real ( kind = dp ), parameter :: p0 =  0.577215664901533D+00
    real ( kind = dp ), parameter :: p1 =  0.844203922187225D+00
    real ( kind = dp ), parameter :: p2 = -0.168860593646662D+00
    real ( kind = dp ), parameter :: p3 = -0.780427615533591D+00
    real ( kind = dp ), parameter :: p4 = -0.402055799310489D+00
    real ( kind = dp ), parameter :: p5 = -0.673562214325671D-01
    real ( kind = dp ), parameter :: p6 = -0.271935708322958D-02
    real ( kind = dp ), parameter :: q1 =  0.288743195473681D+01
    real ( kind = dp ), parameter :: q2 =  0.312755088914843D+01
    real ( kind = dp ), parameter :: q3 =  0.156875193295039D+01
    real ( kind = dp ), parameter :: q4 =  0.361951990101499D+00
    real ( kind = dp ), parameter :: q5 =  0.325038868253937D-01
    real ( kind = dp ), parameter :: q6 =  0.667465618796164D-03
    real ( kind = dp ), parameter :: r0 = 0.422784335098467D+00
    real ( kind = dp ), parameter :: r1 = 0.848044614534529D+00
    real ( kind = dp ), parameter :: r2 = 0.565221050691933D+00
    real ( kind = dp ), parameter :: r3 = 0.156513060486551D+00
    real ( kind = dp ), parameter :: r4 = 0.170502484022650D-01
    real ( kind = dp ), parameter :: r5 = 0.497958207639485D-03
    real ( kind = dp ), parameter :: s1 = 0.124313399877507D+01
    real ( kind = dp ), parameter :: s2 = 0.548042109832463D+00
    real ( kind = dp ), parameter :: s3 = 0.101552187439830D+00
    real ( kind = dp ), parameter :: s4 = 0.713309612391000D-02
    real ( kind = dp ), parameter :: s5 = 0.116165475989616D-03

    if ( a < 0.6D+00 ) then

        top = (((((  &
            p6   &
            * a + p5 ) &
            * a + p4 ) &
            * a + p3 ) &
            * a + p2 ) &
            * a + p1 ) &
            * a + p0

        bot = (((((  &
            q6   &
            * a + q5 ) &
            * a + q4 ) &
            * a + q3 ) &
            * a + q2 ) &
            * a + q1 ) &
            * a + 1.0D+00

        gamma_ln1 = -a * ( top / bot )

    else

        x = ( a - 0.5D+00 ) - 0.5D+00

        top = ((((( r5 * x + r4 ) * x + r3 ) * x + r2 ) * x + r1 ) * x + r0 )

        bot = ((((( s5 * x + s4 ) * x + s3 ) * x + s2 ) * x + s1 ) * x + 1.0D+00 )

        gamma_ln1 = x * ( top / bot )

    end if

    return
    end function gamma_ln1

    function gamma_log ( a )

    !*****************************************************************************80
    !
    !! GAMMA_LOG evaluates ln ( Gamma ( A ) ) for positive A.
    !
    !  Author:
    !
    !    Alfred Morris,
    !    Naval Surface Weapons Center,
    !    Dahlgren, Virginia.
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) A, the argument of the function.
    !    A should be positive.
    !
    !    Output, real ( kind = dp ), GAMMA_LOG, the value of ln ( Gamma ( A ) ).
    !
    implicit none

    real ( kind = dp ) :: a
    real ( kind = dp ) :: gamma_log

    real ( kind = dp ) :: t,w
    integer ( kind = int32 ) :: i,n

    real ( kind = dp ), parameter :: c0 =  0.833333333333333D-01
    real ( kind = dp ), parameter :: c1 = -0.277777777760991D-02
    real ( kind = dp ), parameter :: c2 =  0.793650666825390D-03
    real ( kind = dp ), parameter :: c3 = -0.595202931351870D-03
    real ( kind = dp ), parameter :: c4 =  0.837308034031215D-03
    real ( kind = dp ), parameter :: c5 = -0.165322962780713D-02
    real ( kind = dp ), parameter :: d  =  0.418938533204673D+00

    if ( a <= 0.8D+00 ) then

        gamma_log = gamma_ln1 ( a ) - log ( a )

    else if ( a <= 2.25D+00 ) then

        t = ( a - 0.5D+00 ) - 0.5D+00
        gamma_log = gamma_ln1 ( t )

    else if ( a < 10.0D+00 ) then

        n = a - 1.25D+00
        t = a
        w = 1.0D+00
        do i = 1, n
            t = t - 1.0D+00
            w = t * w
        end do

        gamma_log = gamma_ln1 ( t - 1.0D+00 ) + log ( w )

    else

        t = ( 1.0D+00 / a )**2

        w = ((((( c5 * t + c4 ) * t + c3 ) * t + c2 ) * t + c1 ) * t + c0 ) / a

        gamma_log = ( d + w ) + ( a - 0.5D+00 ) &
            * ( log ( a ) - 1.0D+00 )

    end if

    return
    end function gamma_log

    subroutine gamma_rat1 ( a, x, r, p, q, eps )

    !*****************************************************************************80
    !
    !! GAMMA_RAT1 evaluates the incomplete gamma ratio functions P(A,X) and Q(A,X).
    !
    !  Author:
    !
    !    Armido DiDinato, Alfred Morris
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) A, X, the parameters of the functions.
    !    It is assumed that A <= 1.
    !
    !    Input, real ( kind = dp ) R, the value exp(-X) * X**A / Gamma(A).
    !
    !    Output, real ( kind = dp ) P, Q, the values of P(A,X) and Q(A,X).
    !
    !    Input, real ( kind = dp ) EPS, the tolerance.
    !
    implicit none

    real ( kind = dp ) :: a,x,r,p,q,eps

    real ( kind = dp ) :: a2n,a2nm1,am0,an,an0,b2n,b2nm1,c,cma
    real ( kind = dp ) :: g,h,j,l,sum1,t,tol,w,z

    if ( a * x == 0.0D+00 ) then

        if ( x <= a ) then
            p = 0.0D+00
            q = 1.0D+00
        else
            p = 1.0D+00
            q = 0.0D+00
        end if

        return
    end if

    if ( a == 0.5D+00 ) then

        if ( x < 0.25D+00 ) then
            p = error_f ( sqrt ( x ) )
            q = 0.5D+00 + ( 0.5D+00 - p )
        else
            q = error_fc ( 0, sqrt ( x ) )
            p = 0.5D+00 + ( 0.5D+00 - q )
        end if

        return

    end if
    !
    !  Taylor series for P(A,X)/X**A
    !
    if ( x < 1.1D+00 ) then

        an = 3.0
        c = x
        sum1 = x / ( a + 3.0D+00 )
        tol = 0.1D+00 * eps / ( a + 1.0D+00 )

        do

            an = an + 1.0D+00
            c = -c * ( x / an )
            t = c / ( a + an )
            sum1 = sum1 + t

            if ( abs ( t ) <= tol ) then
                exit
            end if

        end do

        j = a * x * ( ( sum1 / 6.0D+00 - 0.5D+00 &
            / ( a +  2.0D+00  ) ) &
            * x + 1.0D+00 / ( a + 1.0D+00 ) )

        z = a * log ( x )
        h = gam1 ( a )
        g = 1.0D+00 + h

        if ( x < 0.25D+00 ) then
            go to 30
        end if

        if ( a < x / 2.59D+00 ) then
            go to 50
        else
            go to 40
        end if

30      continue

        if ( -0.13394D+00 < z ) then
            go to 50
        end if

40      continue

        w = exp ( z )
        p = w * g * ( 0.5D+00 + ( 0.5D+00 - j ))
        q = 0.5D+00 + ( 0.5D+00 - p )
        return

50      continue

        l = rexp ( z )
        w = 0.5D+00 + ( 0.5D+00 + l )
        q = ( w * j - l ) * g - h

        if  ( q < 0.0D+00 ) then
            p = 1.0D+00
            q = 0.0D+00
        else
            p = 0.5D+00 + ( 0.5D+00 - q )
        end if
        !
        !  Continued fraction expansion.
        !
    else

        a2nm1 = 1.0D+00
        a2n = 1.0D+00
        b2nm1 = x
        b2n = x + ( 1.0D+00 - a )
        c = 1.0D+00

        do

            a2nm1 = x * a2n + c * a2nm1
            b2nm1 = x * b2n + c * b2nm1
            am0 = a2nm1 / b2nm1
            c = c + 1.0D+00
            cma = c - a
            a2n = a2nm1 + cma * a2n
            b2n = b2nm1 + cma * b2n
            an0 = a2n / b2n

            if ( abs ( an0 - am0 ) < eps * an0 ) then
                exit
            end if

        end do

        q = r * an0
        p = 0.5D+00 + ( 0.5D+00 - q )

    end if

    return
    end subroutine gamma_rat1

    subroutine gamma_values ( n_data, x, fx )

    !*****************************************************************************80
    !
    !! GAMMA_VALUES returns some values of the Gamma function.
    !
    !  Definition:
    !
    !    Gamma(Z) = Integral ( 0 <= T < Infinity) T**(Z-1) exp(-T) dT
    !
    !  Recursion:
    !
    !    Gamma(X+1) = X * Gamma(X)
    !
    !  Restrictions:
    !
    !    0 < X ( a software restriction).
    !
    !  Special values:
    !
    !    GAMMA(0.5) = sqrt(PI)
    !
    !    For N a positive integer, GAMMA(N+1) = N!, the standard factorial.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    17 April 2001
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Milton Abramowitz, Irene Stegun,
    !    Handbook of Mathematical Functions,
    !    US Department of Commerce, 1964.
    !
    !  Parameters:
    !
    !    Input/output, integer ( kind = int32 ) N_DATA.  The user sets N_DATA to 0
    !    before the first call.  On each call, the routine increments N_DATA by 1,
    !    and returns the corresponding data; when there is no more data, the
    !    output value of N_DATA will be 0 again.
    !
    !    Output, real ( kind = dp ) X, the argument of the function.
    !
    !    Output, real ( kind = dp ) FX, the value of the function.
    !
    implicit none

    integer ( kind = int32 ) :: n_data
    real ( kind = dp ) :: x,fx

    integer ( kind = int32 ), parameter :: n_max = 18

    real ( kind = dp ), save, dimension ( n_max ) :: fx_vec = (/ &
        4.590845D+00,     2.218160D+00,     1.489192D+00,     1.164230D+00, &
        1.0000000000D+00, 0.9513507699D+00, 0.9181687424D+00, 0.8974706963D+00, &
        0.8872638175D+00, 0.8862269255D+00, 0.8935153493D+00, 0.9086387329D+00, &
        0.9313837710D+00, 0.9617658319D+00, 1.0000000000D+00, 3.6288000D+05, &
        1.2164510D+17,    8.8417620D+30 /)
    real ( kind = dp ), save, dimension ( n_max ) :: x_vec = (/ &
        0.2D+00,  0.4D+00,  0.6D+00,  0.8D+00, &
        1.0D+00,  1.1D+00,  1.2D+00,  1.3D+00, &
        1.4D+00,  1.5D+00,  1.6D+00,  1.7D+00, &
        1.8D+00,  1.9D+00,  2.0D+00, 10.0D+00, &
        20.0D+00, 30.0D+00 /)

    if ( n_data < 0 ) then
        n_data = 0
    end if

    n_data = n_data + 1

    if ( n_max < n_data ) then
        n_data = 0
        x = 0.0D+00
        fx = 0.0D+00
    else
        x = x_vec(n_data)
        fx = fx_vec(n_data)
    end if

    return
    end subroutine gamma_values

    function gsumln ( a, b )

    !*****************************************************************************80
    !
    !! GSUMLN evaluates the function ln(Gamma(A + B)).
    !
    !  Discussion:
    !
    !    GSUMLN is used for 1 <= A <= 2 and 1 <= B <= 2
    !
    !  Author:
    !
    !    Armido DiDinato, Alfred Morris
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) A, B, values whose sum is the argument of
    !    the Gamma function.
    !
    !    Output, real ( kind = dp ) GSUMLN, the value of ln(Gamma(A+B)).
    !
    implicit none

    real ( kind = dp ) :: a,b
    real ( kind = dp ) :: gsumln

    real ( kind = dp ) :: x

    x = a + b - 2.0D+00

    if ( x <= 0.25D+00 ) then
        gsumln = gamma_ln1 ( 1.0D+00 + x )
    else if ( x <= 1.25D+00 ) then
        gsumln = gamma_ln1 ( x ) + alnrel ( x )
    else
        gsumln = gamma_ln1 ( x - 1.0D+00 ) + log ( x * ( 1.0D+00 + x ) )
    end if

    return
    end function gsumln

    function ipmpar ( i )

    !*****************************************************************************80
    !
    !! IPMPAR returns integer machine constants.
    !
    !  Discussion:
    !
    !    Input arguments 1 through 3 are queries about integer arithmetic.
    !    We assume integers are represented in the N-digit, base A form
    !
    !      sign * ( X(N-1)*A^(N-1) + ... + X(1)*A + X(0) )
    !
    !    where 0 <= X(0:N-1) < A.
    !
    !    Then:
    !
    !      IPMPAR(1) = A, the base of integer arithmetic;
    !      IPMPAR(2) = N, the number of base A digits;
    !      IPMPAR(3) = A^N - 1, the largest magnitude.
    !
    !    It is assumed that the single and real ( kind = dp ) floating
    !    point arithmetics have the same base, say B, and that the
    !    nonzero numbers are represented in the form
    !
    !      sign * (B^E) * (X(1)/B + ... + X(M)/B^M)
    !
    !    where X(1:M) is one of { 0, 1,..., B-1 }, and 1 <= X(1) and
    !    EMIN <= E <= EMAX.
    !
    !    Input argument 4 is a query about the base of real arithmetic:
    !
    !      IPMPAR(4) = B, the base of single and real ( kind = dp ) arithmetic.
    !
    !    Input arguments 5 through 7 are queries about single precision
    !    floating point arithmetic:
    !
    !     IPMPAR(5) = M, the number of base B digits for single precision.
    !     IPMPAR(6) = EMIN, the smallest exponent E for single precision.
    !     IPMPAR(7) = EMAX, the largest exponent E for single precision.
    !
    !    Input arguments 8 through 10 are queries about real ( kind = dp )
    !    floating point arithmetic:
    !
    !     IPMPAR(8) = M, the number of base B digits for real ( kind = dp ).
    !     IPMPAR(9) = EMIN, the smallest exponent E for real ( kind = dp ).
    !     IPMPAR(10) = EMAX, the largest exponent E for real ( kind = dp ).
    !
    !  Reference:
    !
    !    Phyllis Fox, Andrew Hall, Norman Schryer,
    !    Algorithm 528:
    !    Framework for a Portable FORTRAN Subroutine Library,
    !    ACM Transactions on Mathematical Software,
    !    Volume 4, 1978, pages 176-188.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = int32 ) I, the index of the desired constant.
    !
    !    Output, integer ( kind = int32 ) IPMPAR, the value of the desired constant.
    !
    implicit none

    integer ( kind = int32 ) :: i
    integer ( kind = int32 ) :: ipmpar

    integer ( kind = int32 ) :: imach(10)

    !
    !     MACHINE CONSTANTS FOR AMDAHL MACHINES.
    !
    !     data imach( 1) /   2 /
    !     data imach( 2) /  31 /
    !     data imach( 3) / 2147483647 /
    !     data imach( 4) /  16 /
    !     data imach( 5) /   6 /
    !     data imach( 6) / -64 /
    !     data imach( 7) /  63 /
    !     data imach( 8) /  14 /
    !     data imach( 9) / -64 /
    !     data imach(10) /  63 /
    !
    !     Machine constants for the AT&T 3B SERIES, AT&T
    !     PC 7300, AND AT&T 6300.
    !
    !     data imach( 1) /     2 /
    !     data imach( 2) /    31 /
    !     data imach( 3) / 2147483647 /
    !     data imach( 4) /     2 /
    !     data imach( 5) /    24 /
    !     data imach( 6) /  -125 /
    !     data imach( 7) /   128 /
    !     data imach( 8) /    53 /
    !     data imach( 9) / -1021 /
    !     data imach(10) /  1024 /
    !
    !     Machine constants for the BURROUGHS 1700 SYSTEM.
    !
    !     data imach( 1) /    2 /
    !     data imach( 2) /   33 /
    !     data imach( 3) / 8589934591 /
    !     data imach( 4) /    2 /
    !     data imach( 5) /   24 /
    !     data imach( 6) / -256 /
    !     data imach( 7) /  255 /
    !     data imach( 8) /   60 /
    !     data imach( 9) / -256 /
    !     data imach(10) /  255 /
    !
    !     Machine constants for the BURROUGHS 5700 SYSTEM.
    !
    !     data imach( 1) /    2 /
    !     data imach( 2) /   39 /
    !     data imach( 3) / 549755813887 /
    !     data imach( 4) /    8 /
    !     data imach( 5) /   13 /
    !     data imach( 6) /  -50 /
    !     data imach( 7) /   76 /
    !     data imach( 8) /   26 /
    !     data imach( 9) /  -50 /
    !     data imach(10) /   76 /
    !
    !     Machine constants for the BURROUGHS 6700/7700 SYSTEMS.
    !
    !     data imach( 1) /      2 /
    !     data imach( 2) /     39 /
    !     data imach( 3) / 549755813887 /
    !     data imach( 4) /      8 /
    !     data imach( 5) /     13 /
    !     data imach( 6) /    -50 /
    !     data imach( 7) /     76 /
    !     data imach( 8) /     26 /
    !     data imach( 9) / -32754 /
    !     data imach(10) /  32780 /
    !
    !     Machine constants for the CDC 6000/7000 SERIES
    !     60 BIT ARITHMETIC, AND THE CDC CYBER 995 64 BIT
    !     ARITHMETIC (NOS OPERATING SYSTEM).
    !
    !     data imach( 1) /    2 /
    !     data imach( 2) /   48 /
    !     data imach( 3) / 281474976710655 /
    !     data imach( 4) /    2 /
    !     data imach( 5) /   48 /
    !     data imach( 6) / -974 /
    !     data imach( 7) / 1070 /
    !     data imach( 8) /   95 /
    !     data imach( 9) / -926 /
    !     data imach(10) / 1070 /
    !
    !     Machine constants for the CDC CYBER 995 64 BIT
    !     ARITHMETIC (NOS/VE OPERATING SYSTEM).
    !
    !     data imach( 1) /     2 /
    !     data imach( 2) /    63 /
    !     data imach( 3) / 9223372036854775807 /
    !     data imach( 4) /     2 /
    !     data imach( 5) /    48 /
    !     data imach( 6) / -4096 /
    !     data imach( 7) /  4095 /
    !     data imach( 8) /    96 /
    !     data imach( 9) / -4096 /
    !     data imach(10) /  4095 /
    !
    !     Machine constants for the CRAY 1, XMP, 2, AND 3.
    !
    !     data imach( 1) /     2 /
    !     data imach( 2) /    63 /
    !     data imach( 3) / 9223372036854775807 /
    !     data imach( 4) /     2 /
    !     data imach( 5) /    47 /
    !     data imach( 6) / -8189 /
    !     data imach( 7) /  8190 /
    !     data imach( 8) /    94 /
    !     data imach( 9) / -8099 /
    !     data imach(10) /  8190 /
    !
    !     Machine constants for the data GENERAL ECLIPSE S/200.
    !
    !     data imach( 1) /    2 /
    !     data imach( 2) /   15 /
    !     data imach( 3) / 32767 /
    !     data imach( 4) /   16 /
    !     data imach( 5) /    6 /
    !     data imach( 6) /  -64 /
    !     data imach( 7) /   63 /
    !     data imach( 8) /   14 /
    !     data imach( 9) /  -64 /
    !     data imach(10) /   63 /
    !
    !     Machine constants for the HARRIS 220.
    !
    !     data imach( 1) /    2 /
    !     data imach( 2) /   23 /
    !     data imach( 3) / 8388607 /
    !     data imach( 4) /    2 /
    !     data imach( 5) /   23 /
    !     data imach( 6) / -127 /
    !     data imach( 7) /  127 /
    !     data imach( 8) /   38 /
    !     data imach( 9) / -127 /
    !     data imach(10) /  127 /
    !
    !     Machine constants for the HONEYWELL 600/6000
    !     AND DPS 8/70 SERIES.
    !
    !     data imach( 1) /    2 /
    !     data imach( 2) /   35 /
    !     data imach( 3) / 34359738367 /
    !     data imach( 4) /    2 /
    !     data imach( 5) /   27 /
    !     data imach( 6) / -127 /
    !     data imach( 7) /  127 /
    !     data imach( 8) /   63 /
    !     data imach( 9) / -127 /
    !     data imach(10) /  127 /
    !
    !     Machine constants for the HP 2100
    !     3 WORD real ( kind = dp ) OPTION WITH FTN4
    !
    !     data imach( 1) /    2 /
    !     data imach( 2) /   15 /
    !     data imach( 3) / 32767 /
    !     data imach( 4) /    2 /
    !     data imach( 5) /   23 /
    !     data imach( 6) / -128 /
    !     data imach( 7) /  127 /
    !     data imach( 8) /   39 /
    !     data imach( 9) / -128 /
    !     data imach(10) /  127 /
    !
    !     Machine constants for the HP 2100
    !     4 WORD real ( kind = dp ) OPTION WITH FTN4
    !
    !     data imach( 1) /    2 /
    !     data imach( 2) /   15 /
    !     data imach( 3) / 32767 /
    !     data imach( 4) /    2 /
    !     data imach( 5) /   23 /
    !     data imach( 6) / -128 /
    !     data imach( 7) /  127 /
    !     data imach( 8) /   55 /
    !     data imach( 9) / -128 /
    !     data imach(10) /  127 /
    !
    !     Machine constants for the HP 9000.
    !
    !     data imach( 1) /     2 /
    !     data imach( 2) /    31 /
    !     data imach( 3) / 2147483647 /
    !     data imach( 4) /     2 /
    !     data imach( 5) /    24 /
    !     data imach( 6) /  -126 /
    !     data imach( 7) /   128 /
    !     data imach( 8) /    53 /
    !     data imach( 9) / -1021 /
    !     data imach(10) /  1024 /
    !
    !     Machine constants for the IBM 360/370 SERIES,
    !     THE ICL 2900, THE ITEL AS/6, THE XEROX SIGMA
    !     5/7/9 AND THE SEL SYSTEMS 85/86.
    !
    !     data imach( 1) /    2 /
    !     data imach( 2) /   31 /
    !     data imach( 3) / 2147483647 /
    !     data imach( 4) /   16 /
    !     data imach( 5) /    6 /
    !     data imach( 6) /  -64 /
    !     data imach( 7) /   63 /
    !     data imach( 8) /   14 /
    !     data imach( 9) /  -64 /
    !     data imach(10) /   63 /
    !
    !     Machine constants for the IBM PC.
    !
    !      data imach(1)/2/
    !      data imach(2)/31/
    !      data imach(3)/2147483647/
    !      data imach(4)/2/
    !      data imach(5)/24/
    !      data imach(6)/-125/
    !      data imach(7)/128/
    !      data imach(8)/53/
    !      data imach(9)/-1021/
    !      data imach(10)/1024/
    !
    !     Machine constants for the MACINTOSH II - ABSOFT
    !     MACFORTRAN II.
    !
    !     data imach( 1) /     2 /
    !     data imach( 2) /    31 /
    !     data imach( 3) / 2147483647 /
    !     data imach( 4) /     2 /
    !     data imach( 5) /    24 /
    !     data imach( 6) /  -125 /
    !     data imach( 7) /   128 /
    !     data imach( 8) /    53 /
    !     data imach( 9) / -1021 /
    !     data imach(10) /  1024 /
    !
    !     Machine constants for the MICROVAX - VMS FORTRAN.
    !
    !     data imach( 1) /    2 /
    !     data imach( 2) /   31 /
    !     data imach( 3) / 2147483647 /
    !     data imach( 4) /    2 /
    !     data imach( 5) /   24 /
    !     data imach( 6) / -127 /
    !     data imach( 7) /  127 /
    !     data imach( 8) /   56 /
    !     data imach( 9) / -127 /
    !     data imach(10) /  127 /
    !
    !     Machine constants for the PDP-11 FORTRAN SUPPORTING
    !     32-BIT integer ARITHMETIC.
    !
    !     data imach( 1) /    2 /
    !     data imach( 2) /   31 /
    !     data imach( 3) / 2147483647 /
    !     data imach( 4) /    2 /
    !     data imach( 5) /   24 /
    !     data imach( 6) / -127 /
    !     data imach( 7) /  127 /
    !     data imach( 8) /   56 /
    !     data imach( 9) / -127 /
    !     data imach(10) /  127 /
    !
    !     Machine constants for the SEQUENT BALANCE 8000.
    !
    !     data imach( 1) /     2 /
    !     data imach( 2) /    31 /
    !     data imach( 3) / 2147483647 /
    !     data imach( 4) /     2 /
    !     data imach( 5) /    24 /
    !     data imach( 6) /  -125 /
    !     data imach( 7) /   128 /
    !     data imach( 8) /    53 /
    !     data imach( 9) / -1021 /
    !     data imach(10) /  1024 /
    !
    !     Machine constants for the SILICON GRAPHICS IRIS-4D
    !     SERIES (MIPS R3000 PROCESSOR).
    !
    !     data imach( 1) /     2 /
    !     data imach( 2) /    31 /
    !     data imach( 3) / 2147483647 /
    !     data imach( 4) /     2 /
    !     data imach( 5) /    24 /
    !     data imach( 6) /  -125 /
    !     data imach( 7) /   128 /
    !     data imach( 8) /    53 /
    !     data imach( 9) / -1021 /
    !     data imach(10) /  1024 /
    !
    !     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
    !     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
    !     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
    !
    data imach( 1) /     2 /
    data imach( 2) /    31 /
    data imach( 3) / 2147483647 /
    data imach( 4) /     2 /
    data imach( 5) /    24 /
    data imach( 6) /  -125 /
    data imach( 7) /   128 /
    data imach( 8) /    53 /
    data imach( 9) / -1021 /
    data imach(10) /  1024 /
    !
    !     Machine constants for the UNIVAC 1100 SERIES.
    !
    !     data imach( 1) /    2 /
    !     data imach( 2) /   35 /
    !     data imach( 3) / 34359738367 /
    !     data imach( 4) /    2 /
    !     data imach( 5) /   27 /
    !     data imach( 6) / -128 /
    !     data imach( 7) /  127 /
    !     data imach( 8) /   60 /
    !     data imach( 9) /-1024 /
    !     data imach(10) / 1023 /
    !
    !     Machine constants for the VAX 11/780.
    !
    !     data imach( 1) /    2 /
    !     data imach( 2) /   31 /
    !     data imach( 3) / 2147483647 /
    !     data imach( 4) /    2 /
    !     data imach( 5) /   24 /
    !     data imach( 6) / -127 /
    !     data imach( 7) /  127 /
    !     data imach( 8) /   56 /
    !     data imach( 9) / -127 /
    !     data imach(10) /  127 /
    !

    ! alternatively, use Fortran intrinsics
    imach(1)  = radix(1)           ! A = base of integer arithmetic
    imach(2)  = digits(1)          ! # of base A digits for integers
    imach(3)  = huge(1)            ! A^(N-1) : maximum integer value
    imach(4)  = radix(2.0d0)       ! B = base of real arithmetic
    imach(5)  = digits(1.0)        ! number of base B digits for single precision
    imach(6)  = minexponent(1.0)   ! min exponent for single precision
    imach(7)  = maxexponent(1.0)   ! max exponent for single precision
    imach(8)  = digits(1.0d0)      ! number of base B digits for double precision
    imach(9)  = minexponent(1.0d0) ! min exponent for double precision
    imach(10) = maxexponent(1.0d0) ! max exponent for double precision

    ipmpar = imach(i)

    return
    end function ipmpar

    subroutine negative_binomial_cdf_values ( n_data, f, s, p, cdf )

    !*****************************************************************************80
    !
    !! NEGATIVE_BINOMIAL_CDF_VALUES returns values of the negative binomial CDF.
    !
    !  Discussion:
    !
    !    Assume that a coin has a probability P of coming up heads on
    !    any one trial.  Suppose that we plan to flip the coin until we
    !    achieve a total of S heads.  If we let F represent the number of
    !    tails that occur in this process, then the value of F satisfies
    !    a negative binomial PDF:
    !
    !      PDF(F,S,P) = Choose ( F from F+S-1 ) * P**S * (1-P)**F
    !
    !    The negative binomial CDF is the probability that there are F or
    !    fewer failures upon the attainment of the S-th success.  Thus,
    !
    !      CDF(F,S,P) = sum ( 0 <= G <= F ) PDF(G,S,P)
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    07 June 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    FC Powell,
    !    Statistical Tables for Sociology, Biology and Physical Sciences,
    !    Cambridge University Press, 1982.
    !
    !  Parameters:
    !
    !    Input/output, integer ( kind = int32 ) N_DATA.  The user sets N_DATA to 0
    !    before the first call.  On each call, the routine increments N_DATA by 1,
    !    and returns the corresponding data; when there is no more data, the
    !    output value of N_DATA will be 0 again.
    !
    !    Output, integer ( kind = int32 ) F, the maximum number of failures.
    !
    !    Output, integer ( kind = int32 ) S, the number of successes.
    !
    !    Output, real ( kind = dp ) P, the probability of a success on one trial.
    !
    !    Output, real ( kind = dp ) CDF, the probability of at most F failures
    !    before the S-th success.
    !
    implicit none

    integer ( kind = int32 ) :: n_data,f,s
    real ( kind = dp ) :: p,cdf

    integer ( kind = int32 ), parameter :: n_max = 27

    real ( kind = dp ), save, dimension ( n_max ) :: cdf_vec = (/ &
        0.6367D+00, 0.3633D+00, 0.1445D+00, &
        0.5000D+00, 0.2266D+00, 0.0625D+00, &
        0.3438D+00, 0.1094D+00, 0.0156D+00, &
        0.1792D+00, 0.0410D+00, 0.0041D+00, &
        0.0705D+00, 0.0109D+00, 0.0007D+00, &
        0.9862D+00, 0.9150D+00, 0.7472D+00, &
        0.8499D+00, 0.5497D+00, 0.2662D+00, &
        0.6513D+00, 0.2639D+00, 0.0702D+00, &
        1.0000D+00, 0.0199D+00, 0.0001D+00 /)
    integer ( kind = int32 ), save, dimension ( n_max ) :: f_vec = (/ &
        4,  3,  2, &
        3,  2,  1, &
        2,  1,  0, &
        2,  1,  0, &
        2,  1,  0, &
        11, 10,  9, &
        17, 16, 15, &
        9,  8,  7, &
        2,  1,  0 /)
    real ( kind = dp ), save, dimension ( n_max ) :: p_vec = (/ &
        0.50D+00, 0.50D+00, 0.50D+00, &
        0.50D+00, 0.50D+00, 0.50D+00, &
        0.50D+00, 0.50D+00, 0.50D+00, &
        0.40D+00, 0.40D+00, 0.40D+00, &
        0.30D+00, 0.30D+00, 0.30D+00, &
        0.30D+00, 0.30D+00, 0.30D+00, &
        0.10D+00, 0.10D+00, 0.10D+00, &
        0.10D+00, 0.10D+00, 0.10D+00, &
        0.01D+00, 0.01D+00, 0.01D+00 /)
    integer ( kind = int32 ), save, dimension ( n_max ) :: s_vec = (/ &
        4, 5, 6, &
        4, 5, 6, &
        4, 5, 6, &
        4, 5, 6, &
        4, 5, 6, &
        1, 2, 3, &
        1, 2, 3, &
        1, 2, 3, &
        0, 1, 2 /)

    if ( n_data < 0 ) then
        n_data = 0
    end if

    n_data = n_data + 1

    if ( n_max < n_data ) then
        n_data = 0
        f = 0
        s = 0
        p = 0.0D+00
        cdf = 0.0D+00
    else
        f = f_vec(n_data)
        s = s_vec(n_data)
        p = p_vec(n_data)
        cdf = cdf_vec(n_data)
    end if

    return
    end subroutine negative_binomial_cdf_values

    subroutine normal_01_cdf_values ( n_data, x, fx )

    !*****************************************************************************80
    !
    !! NORMAL_01_CDF_VALUES returns some values of the Normal 01 CDF.
    !
    !  Discussion:
    !
    !    In Mathematica, the function can be evaluated by:
    !
    !      Needs["Statistics`ContinuousDistributions`"]
    !      dist = NormalDistribution [ 0, 1 ]
    !      CDF [ dist, x ]
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    28 August 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Milton Abramowitz, Irene Stegun,
    !    Handbook of Mathematical Functions,
    !    US Department of Commerce, 1964.
    !
    !    Stephen Wolfram,
    !    The Mathematica Book,
    !    Fourth Edition,
    !    Wolfram Media / Cambridge University Press, 1999.
    !
    !  Parameters:
    !
    !    Input/output, integer ( kind = int32 ) N_DATA.  The user sets N_DATA to 0
    !    before the first call.  On each call, the routine increments N_DATA by 1,
    !    and returns the corresponding data; when there is no more data, the
    !    output value of N_DATA will be 0 again.
    !
    !    Output, real ( kind = dp ) X, the argument of the function.
    !
    !    Output, real ( kind = dp ) FX, the value of the function.
    !
    implicit none

    integer ( kind = int32 ) :: n_data
    real ( kind = dp ) :: x,fx

    integer ( kind = int32 ), parameter :: n_max = 17

    real ( kind = dp ), save, dimension ( n_max ) :: fx_vec = (/ &
        0.5000000000000000D+00, &
        0.5398278372770290D+00, &
        0.5792597094391030D+00, &
        0.6179114221889526D+00, &
        0.6554217416103242D+00, &
        0.6914624612740131D+00, &
        0.7257468822499270D+00, &
        0.7580363477769270D+00, &
        0.7881446014166033D+00, &
        0.8159398746532405D+00, &
        0.8413447460685429D+00, &
        0.9331927987311419D+00, &
        0.9772498680518208D+00, &
        0.9937903346742239D+00, &
        0.9986501019683699D+00, &
        0.9997673709209645D+00, &
        0.9999683287581669D+00 /)
    real ( kind = dp ), save, dimension ( n_max ) :: x_vec = (/ &
        0.0000000000000000D+00, &
        0.1000000000000000D+00, &
        0.2000000000000000D+00, &
        0.3000000000000000D+00, &
        0.4000000000000000D+00, &
        0.5000000000000000D+00, &
        0.6000000000000000D+00, &
        0.7000000000000000D+00, &
        0.8000000000000000D+00, &
        0.9000000000000000D+00, &
        0.1000000000000000D+01, &
        0.1500000000000000D+01, &
        0.2000000000000000D+01, &
        0.2500000000000000D+01, &
        0.3000000000000000D+01, &
        0.3500000000000000D+01, &
        0.4000000000000000D+01 /)

    if ( n_data < 0 ) then
        n_data = 0
    end if

    n_data = n_data + 1

    if ( n_max < n_data ) then
        n_data = 0
        x = 0.0D+00
        fx = 0.0D+00
    else
        x = x_vec(n_data)
        fx = fx_vec(n_data)
    end if

    return
    end subroutine normal_01_cdf_values

    subroutine normal_cdf_values ( n_data, mu, sigma, x, fx )

    !*****************************************************************************80
    !
    !! NORMAL_CDF_VALUES returns some values of the Normal CDF.
    !
    !  Discussion:
    !
    !    In Mathematica, the function can be evaluated by:
    !
    !      Needs["Statistics`ContinuousDistributions`"]
    !      dist = NormalDistribution [ mu, sigma ]
    !      CDF [ dist, x ]
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 August 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Milton Abramowitz, Irene Stegun,
    !    Handbook of Mathematical Functions,
    !    US Department of Commerce, 1964.
    !
    !    Stephen Wolfram,
    !    The Mathematica Book,
    !    Fourth Edition,
    !    Wolfram Media / Cambridge University Press, 1999.
    !
    !  Parameters:
    !
    !    Input/output, integer ( kind = int32 ) N_DATA.  The user sets N_DATA to 0
    !    before the first call.  On each call, the routine increments N_DATA by 1,
    !    and returns the corresponding data; when there is no more data, the
    !    output value of N_DATA will be 0 again.
    !
    !    Output, real ( kind = dp ) MU, the mean of the distribution.
    !
    !    Output, real ( kind = dp ) SIGMA, the variance of the distribution.
    !
    !    Output, real ( kind = dp ) X, the argument of the function.
    !
    !    Output, real ( kind = dp ) FX, the value of the function.
    !
    implicit none

    integer ( kind = int32 ) :: n_data
    real ( kind = dp ) :: mu,sigma,x,fx

    integer ( kind = int32 ), parameter :: n_max = 12

    real ( kind = dp ), save, dimension ( n_max ) :: fx_vec = (/ &
        0.5000000000000000D+00, &
        0.9772498680518208D+00, &
        0.9999683287581669D+00, &
        0.9999999990134124D+00, &
        0.6914624612740131D+00, &
        0.6305586598182364D+00, &
        0.5987063256829237D+00, &
        0.5792597094391030D+00, &
        0.6914624612740131D+00, &
        0.5000000000000000D+00, &
        0.3085375387259869D+00, &
        0.1586552539314571D+00 /)
    real ( kind = dp ), save, dimension ( n_max ) :: mu_vec = (/ &
        0.1000000000000000D+01, &
        0.1000000000000000D+01, &
        0.1000000000000000D+01, &
        0.1000000000000000D+01, &
        0.1000000000000000D+01, &
        0.1000000000000000D+01, &
        0.1000000000000000D+01, &
        0.1000000000000000D+01, &
        0.2000000000000000D+01, &
        0.3000000000000000D+01, &
        0.4000000000000000D+01, &
        0.5000000000000000D+01 /)
    real ( kind = dp ), save, dimension ( n_max ) :: sigma_vec = (/ &
        0.5000000000000000D+00, &
        0.5000000000000000D+00, &
        0.5000000000000000D+00, &
        0.5000000000000000D+00, &
        0.2000000000000000D+01, &
        0.3000000000000000D+01, &
        0.4000000000000000D+01, &
        0.5000000000000000D+01, &
        0.2000000000000000D+01, &
        0.2000000000000000D+01, &
        0.2000000000000000D+01, &
        0.2000000000000000D+01 /)
    real ( kind = dp ), save, dimension ( n_max ) :: x_vec = (/ &
        0.1000000000000000D+01, &
        0.2000000000000000D+01, &
        0.3000000000000000D+01, &
        0.4000000000000000D+01, &
        0.2000000000000000D+01, &
        0.2000000000000000D+01, &
        0.2000000000000000D+01, &
        0.2000000000000000D+01, &
        0.3000000000000000D+01, &
        0.3000000000000000D+01, &
        0.3000000000000000D+01, &
        0.3000000000000000D+01 /)

    if ( n_data < 0 ) then
        n_data = 0
    end if

    n_data = n_data + 1

    if ( n_max < n_data ) then
        n_data = 0
        mu = 0.0D+00
        sigma = 0.0D+00
        x = 0.0D+00
        fx = 0.0D+00
    else
        mu = mu_vec(n_data)
        sigma = sigma_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
    end if

    return
    end subroutine normal_cdf_values

    subroutine poisson_cdf_values ( n_data, a, x, fx )

    !*****************************************************************************80
    !
    !! POISSON_CDF_VALUES returns some values of the Poisson CDF.
    !
    !  Discussion:
    !
    !    CDF(X)(A) is the probability of at most X successes in unit time,
    !    given that the expected mean number of successes is A.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    28 May 2001
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Milton Abramowitz, Irene Stegun,
    !    Handbook of Mathematical Functions,
    !    US Department of Commerce, 1964.
    !
    !    Daniel Zwillinger,
    !    CRC Standard Mathematical Tables and Formulae,
    !    30th Edition, CRC Press, 1996, pages 653-658.
    !
    !  Parameters:
    !
    !    Input/output, integer ( kind = int32 ) N_DATA.  The user sets N_DATA to 0
    !    before the first call.  On each call, the routine increments N_DATA by 1,
    !    and returns the corresponding data; when there is no more data, the
    !    output value of N_DATA will be 0 again.
    !
    !    Output, real ( kind = dp ) A, integer X, the arguments of the function.
    !
    !    Output, real ( kind = dp ) FX, the value of the function.
    !
    implicit none

    integer ( kind = int32 ) :: n_data
    real ( kind = dp ) :: a
    integer ( kind = int32 ) :: x
    real ( kind = dp ) :: fx

    integer ( kind = int32 ), parameter :: n_max = 21

    real ( kind = dp ), save, dimension ( n_max ) :: a_vec = (/ &
        0.02D+00, 0.10D+00, 0.10D+00, 0.50D+00, &
        0.50D+00, 0.50D+00, 1.00D+00, 1.00D+00, &
        1.00D+00, 1.00D+00, 2.00D+00, 2.00D+00, &
        2.00D+00, 2.00D+00, 5.00D+00, 5.00D+00, &
        5.00D+00, 5.00D+00, 5.00D+00, 5.00D+00, &
        5.00D+00 /)
    real ( kind = dp ), save, dimension ( n_max ) :: fx_vec = (/ &
        0.980D+00, 0.905D+00, 0.995D+00, 0.607D+00, &
        0.910D+00, 0.986D+00, 0.368D+00, 0.736D+00, &
        0.920D+00, 0.981D+00, 0.135D+00, 0.406D+00, &
        0.677D+00, 0.857D+00, 0.007D+00, 0.040D+00, &
        0.125D+00, 0.265D+00, 0.441D+00, 0.616D+00, &
        0.762D+00 /)
    integer ( kind = int32 ), save, dimension ( n_max ) :: x_vec = (/ &
        0, 0, 1, 0, &
        1, 2, 0, 1, &
        2, 3, 0, 1, &
        2, 3, 0, 1, &
        2, 3, 4, 5, &
        6 /)

    if ( n_data < 0 ) then
        n_data = 0
    end if

    n_data = n_data + 1

    if ( n_max < n_data ) then
        n_data = 0
        a = 0.0D+00
        x = 0
        fx = 0.0D+00
    else
        a = a_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
    end if

    return
    end subroutine poisson_cdf_values

    function psi ( xx )

    !*****************************************************************************80
    !
    !! PSI evaluates the psi or digamma function, d/dx ln(gamma(x)).
    !
    !  Discussion:
    !
    !    The main computation involves evaluation of rational Chebyshev
    !    approximations.  PSI was written at Argonne National Laboratory
    !    for FUNPACK, and subsequently modified by A. H. Morris of NSWC.
    !
    !  Reference:
    !
    !    William Cody, Anthony Strecok, Henry Thacher,
    !    Chebyshev Approximations for the Psi Function,
    !    Mathematics of Computation,
    !    Volume 27, 1973, pages 123-127.
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) XX, the argument of the psi function.
    !
    !    Output, real ( kind = dp ) PSI, the value of the psi function.  PSI
    !    is assigned the value 0 when the psi function is undefined.
    !
    implicit none

    real ( kind = dp ) :: xx
    real ( kind = dp ) :: psi

    real ( kind = dp ) :: aug,den,sgn,upper,w,x,xmax1,xmx0,xsmall,z
    integer ( kind = int32 ) :: i,m,n,nq

    real ( kind = dp ), parameter :: dx0 = 1.461632144968362341262659542325721325D+00
    real ( kind = dp ), parameter :: piov4 = 0.785398163397448D+00

    real ( kind = dp ), parameter, dimension ( 7 ) :: p1 = (/ &
        0.895385022981970D-02, &
        0.477762828042627D+01, &
        0.142441585084029D+03, &
        0.118645200713425D+04, &
        0.363351846806499D+04, &
        0.413810161269013D+04, &
        0.130560269827897D+04/)
    real ( kind = dp ), dimension ( 4 ) :: p2 = (/ &
        -0.212940445131011D+01, &
        -0.701677227766759D+01, &
        -0.448616543918019D+01, &
        -0.648157123766197D+00 /)

    !
    !  Coefficients for rational approximation of
    !  PSI(X) / (X - X0),  0.5D+00 <= X <= 3.0D+00
    !
    real ( kind = dp ), dimension ( 6 ) :: q1 = (/ &
        0.448452573429826D+02, &
        0.520752771467162D+03, &
        0.221000799247830D+04, &
        0.364127349079381D+04, &
        0.190831076596300D+04, &
        0.691091682714533D-05 /)
    real ( kind = dp ), dimension ( 4 ) :: q2 = (/ &
        0.322703493791143D+02, &
        0.892920700481861D+02, &
        0.546117738103215D+02, &
        0.777788548522962D+01 /)

    !
    !  XMAX1 is the largest positive floating point constant with entirely
    !  integer representation.  It is also used as negative of lower bound
    !  on acceptable negative arguments and as the positive argument beyond which
    !  psi may be represented as LOG(X).
    !
    xmax1 = real ( ipmpar(3), kind = dp )
    xmax1 = min ( xmax1, 1.0D+00 / epsilon ( xmax1 ) )
    !
    !  XSMALL is the absolute argument below which PI*COTAN(PI*X)
    !  may be represented by 1/X.
    !
    xsmall = 1.0D-09

    x = xx
    aug = 0.0D+00

    if ( x == 0.0D+00 ) then
        psi = 0.0D+00
        return
    end if
    !
    !  X < 0.5,  Use reflection formula PSI(1-X) = PSI(X) + PI * COTAN(PI*X)
    !
    if ( x < 0.5D+00 ) then
        !
        !  0 < ABS ( X ) <= XSMALL.  Use 1/X as a substitute for PI*COTAN(PI*X)
        !
        if ( abs ( x ) <= xsmall ) then
            aug = -1.0D+00 / x
            go to 40
        end if
        !
        !  Reduction of argument for cotangent.
        !
        w = -x
        sgn = piov4

        if ( w <= 0.0D+00 ) then
            w = -w
            sgn = -sgn
        end if
        !
        !  Make an error exit if X <= -XMAX1
        !
        if ( xmax1 <= w ) then
            psi = 0.0D+00
            return
        end if

        nq = int ( w )
        w = w - real ( nq, kind = dp )
        nq = int ( w * 4.0D+00 )
        w = 4.0D+00 * ( w - real ( nq, kind = dp ) * 0.25D+00 )
        !
        !  W is now related to the fractional part of 4.0D+00 * X.
        !  Adjust argument to correspond to values in first
        !  quadrant and determine sign.
        !
        n = nq / 2
        if ( n + n /= nq ) then
            w = 1.0D+00 - w
        end if

        z = piov4 * w
        m = n / 2

        if ( m + m /= n ) then
            sgn = -sgn
        end if
        !
        !  Determine final value for -PI * COTAN(PI*X).
        !
        n = ( nq + 1 ) / 2
        m = n / 2
        m = m + m

        if ( m == n ) then

            if ( z == 0.0D+00 ) then
                psi = 0.0D+00
                return
            end if

            aug = 4.0D+00 * sgn * ( cos(z) / sin(z) )

        else

            aug = 4.0D+00 * sgn * ( sin(z) / cos(z) )

        end if

40      continue

        x = 1.0D+00 - x

    end if
    !
    !  0.5 <= X <= 3
    !
    if ( x <= 3.0D+00 ) then

        den = x
        upper = p1(1) * x

        do i = 1, 5
            den = ( den + q1(i) ) * x
            upper = ( upper + p1(i+1) ) * x
        end do

        den = ( upper + p1(7) ) / ( den + q1(6) )
        xmx0 = real ( x, kind = dp ) - dx0
        psi = den * xmx0 + aug
        !
        !  3 < X < XMAX1
        !
    else if ( x < xmax1 ) then

        w = 1.0D+00 / x**2
        den = w
        upper = p2(1) * w

        do i = 1, 3
            den = ( den + q2(i) ) * w
            upper = ( upper + p2(i+1) ) * w
        end do

        aug = upper / ( den + q2(4) ) - 0.5D+00 / x + aug
        psi = aug + log ( x )
        !
        !  XMAX1 <= X
        !
    else

        psi = aug + log ( x )

    end if

    return
    end function psi

    subroutine psi_values ( n_data, x, fx )

    !*****************************************************************************80
    !
    !! PSI_VALUES returns some values of the Psi or Digamma function.
    !
    !  Discussion:
    !
    !    PSI(X) = d LN ( GAMMA ( X ) ) / d X = GAMMA'(X) / GAMMA(X)
    !
    !    PSI(1) = - Euler's constant.
    !
    !    PSI(X+1) = PSI(X) + 1 / X.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    17 May 2001
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Milton Abramowitz, Irene Stegun,
    !    Handbook of Mathematical Functions,
    !    US Department of Commerce, 1964.
    !
    !  Parameters:
    !
    !    Input/output, integer ( kind = int32 ) N_DATA.  The user sets N_DATA to 0
    !    before the first call.  On each call, the routine increments N_DATA by 1,
    !    and returns the corresponding data; when there is no more data, the
    !    output value of N_DATA will be 0 again.
    !
    !    Output, real ( kind = dp ) X, the argument of the function.
    !
    !    Output, real ( kind = dp ) FX, the value of the function.
    !
    implicit none

    integer ( kind = int32 ) :: n_data
    real ( kind = dp ) :: x,fx

    integer ( kind = int32 ), parameter :: n_max = 11

    real ( kind = dp ), save, dimension ( n_max ) :: fx_vec = (/ &
        -0.5772156649D+00, -0.4237549404D+00, -0.2890398966D+00, &
        -0.1691908889D+00, -0.0613845446D+00, -0.0364899740D+00, &
        0.1260474528D+00,  0.2085478749D+00,  0.2849914333D+00, &
        0.3561841612D+00,  0.4227843351D+00 /)
    real ( kind = dp ), save, dimension ( n_max ) :: x_vec = (/ &
        1.0D+00,  1.1D+00,  1.2D+00,  &
        1.3D+00,  1.4D+00,  1.5D+00,  &
        1.6D+00,  1.7D+00,  1.8D+00,  &
        1.9D+00,  2.0D+00 /)

    if ( n_data < 0 ) then
        n_data = 0
    end if

    n_data = n_data + 1

    if ( n_max < n_data ) then
        n_data = 0
        x = 0.0D+00
        fx = 0.0D+00
    else
        x = x_vec(n_data)
        fx = fx_vec(n_data)
    end if

    return
    end subroutine psi_values

    subroutine r8_swap ( x, y )

    !*****************************************************************************80
    !
    !! R8_SWAP swaps two R8 values.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    01 May 2000
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input/output, real ( kind = dp ) X, Y.  On output, the values of X and
    !    Y have been interchanged.
    !
    implicit none

    real ( kind = dp ) :: x,y

    real ( kind = dp ) :: z

    z = x
    x = y
    y = z

    return
    end subroutine r8_swap

    function rcomp ( a, x )

    !*****************************************************************************80
    !
    !! RCOMP evaluates exp(-X) * X**A / Gamma(A).
    !
    !  Author:
    !
    !    Armido DiDinato, Alfred Morris
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) A, X, arguments of the quantity to be computed.
    !
    !    Output, real ( kind = dp ) RCOMP, the value of exp(-X) * X**A / Gamma(A).
    !
    !  Local parameters:
    !
    !    RT2PIN = 1/SQRT(2*PI)
    !
    implicit none

    real ( kind = dp ) :: a,x
    real ( kind = dp ) :: rcomp

    real ( kind = dp ) :: t,t1,u

    real ( kind = dp ), parameter :: rt2pin = 0.398942280401433D+00


    if ( a < 20.0D+00 ) then

        t = a * log ( x ) - x

        if ( a < 1.0D+00 ) then
            rcomp = ( a * exp ( t ) ) * ( 1.0D+00 + gam1 ( a ) )
        else
            rcomp = exp ( t ) / gamma_cdflib ( a )
        end if

    else

        u = x / a

        if ( u == 0.0D+00 ) then
            rcomp = 0.0D+00
        else
            t = ( 1.0D+00 / a )**2
            t1 = ((( 0.75D+00 * t - 1.0D+00 ) * t + 3.5D+00 ) * t - 105.0D+00 ) &
                / ( a * 1260.0D+00 )
            t1 = t1 - a * rlog ( u )
            rcomp = rt2pin * sqrt ( a ) * exp ( t1 )
        end if

    end if

    return
    end function rcomp

    function rexp ( x )

    !*****************************************************************************80
    !
    !! REXP evaluates the function EXP(X) - 1.
    !
    !  Modified:
    !
    !    09 December 1999
    !
    !  Author:
    !
    !    Armido DiDinato, Alfred Morris
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) X, the argument of the function.
    !
    !    Output, real ( kind = dp ) REXP, the value of EXP(X)-1.
    !
    implicit none

    real ( kind = dp ) :: x
    real ( kind = dp ) :: rexp

    real ( kind = dp ) :: w

    real ( kind = dp ), parameter :: p1 =  0.914041914819518D-09
    real ( kind = dp ), parameter :: p2 =  0.238082361044469D-01
    real ( kind = dp ), parameter :: q1 = -0.499999999085958D+00
    real ( kind = dp ), parameter :: q2 =  0.107141568980644D+00
    real ( kind = dp ), parameter :: q3 = -0.119041179760821D-01
    real ( kind = dp ), parameter :: q4 =  0.595130811860248D-03

    if ( abs ( x ) <= 0.15D+00 ) then

        rexp = x * ( ( ( p2 * x + p1 ) * x + 1.0D+00 ) &
            / ( ( ( ( q4 * x + q3 ) * x + q2 ) * x + q1 ) * x + 1.0D+00 ) )

    else

        w = exp ( x )

        if ( x <= 0.0D+00 ) then
            rexp = ( w - 0.5D+00 ) - 0.5D+00
        else
            rexp = w * ( 0.5D+00 + ( 0.5D+00 - 1.0D+00 / w ) )
        end if

    end if

    return
    end function rexp

    function rlog ( x )

    !*****************************************************************************80
    !
    !! RLOG computes X - 1 - LN(X).
    !
    !  Modified:
    !
    !    06 August 2004
    !
    !  Author:
    !
    !    Armido DiDinato, Alfred Morris
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) X, the argument of the function.
    !
    !    Output, real ( kind = dp ) RLOG, the value of the function.
    !
    implicit none

    real ( kind = dp ) :: x
    real ( kind = dp ) :: rlog

    real ( kind = dp ) :: r,t,u,w,w1

    real ( kind = dp ), parameter :: a  =  0.566749439387324D-01
    real ( kind = dp ), parameter :: b  =  0.456512608815524D-01
    real ( kind = dp ), parameter :: half = 0.5D+00
    real ( kind = dp ), parameter :: p0 =  0.333333333333333D+00
    real ( kind = dp ), parameter :: p1 = -0.224696413112536D+00
    real ( kind = dp ), parameter :: p2 =  0.620886815375787D-02
    real ( kind = dp ), parameter :: q1 = -0.127408923933623D+01
    real ( kind = dp ), parameter :: q2 =  0.354508718369557D+00
    real ( kind = dp ), parameter :: two =  2.0D+00

    if ( x < 0.61D+00 ) then

        r = ( x - 0.5D+00 ) - 0.5D+00
        rlog = r - log ( x )

    else if ( x < 1.57D+00 ) then

        if ( x < 0.82D+00 ) then

            u = x - 0.7D+00
            u = u / 0.7D+00
            w1 = a - u * 0.3D+00

        else if ( x < 1.18D+00 ) then

            u = ( x - half ) - half
            w1 = 0.0D+00

        else if ( x < 1.57D+00 ) then

            u = 0.75D+00 * x - 1.0D+00
            w1 = b + u / 3.0D+00

        end if

        r = u / ( u + two )
        t = r * r
        w = ( ( p2 * t + p1 ) * t + p0 ) / ( ( q2 * t + q1 ) * t + 1.0D+00 )
        rlog = two * t * ( 1.0D+00 / ( 1.0D+00 - r ) - r * w ) + w1

    else if ( 1.57D+00 <= x ) then

        r = ( x - half ) - half
        rlog = r - log ( x )

    end if

    return
    end function rlog

    function rlog1 ( x )

    !*****************************************************************************80
    !
    !! RLOG1 evaluates the function X - ln ( 1 + X ).
    !
    !  Author:
    !
    !    Armido DiDinato, Alfred Morris
    !
    !  Reference:
    !
    !    Armido DiDinato, Alfred Morris,
    !    Algorithm 708:
    !    Significant Digit Computation of the Incomplete Beta Function Ratios,
    !    ACM Transactions on Mathematical Software,
    !    Volume 18, 1993, pages 360-373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) X, the argument.
    !
    !    Output, real ( kind = dp ) RLOG1, the value of X - ln ( 1 + X ).
    !
    implicit none

    real ( kind = dp ) :: x
    real ( kind = dp ) :: rlog1

    real ( kind = dp ) :: h,r,t, w,w1

    real ( kind = dp ), parameter :: a = 0.566749439387324D-01
    real ( kind = dp ), parameter :: b = 0.456512608815524D-01
    real ( kind = dp ), parameter :: p0 = 0.333333333333333D+00
    real ( kind = dp ), parameter :: p1 = -0.224696413112536D+00
    real ( kind = dp ), parameter :: p2 = 0.620886815375787D-02
    real ( kind = dp ), parameter :: q1 = -0.127408923933623D+01
    real ( kind = dp ), parameter :: q2 = 0.354508718369557D+00
    real ( kind = dp ), parameter :: half = 0.5D+00
    real ( kind = dp ), parameter :: two =  2.0D+00

    if ( x < -0.39D+00 ) then

        w = ( x + half ) + half
        rlog1 = x - log ( w )

    else if ( x < -0.18D+00 ) then

        h = x + 0.3D+00
        h = h / 0.7D+00
        w1 = a - h * 0.3D+00

        r = h / ( h + 2.0D+00 )
        t = r * r
        w = ( ( p2 * t + p1 ) * t + p0 ) / ( ( q2 * t + q1 ) * t + 1.0D+00 )
        rlog1 = two * t * ( 1.0D+00 / ( 1.0D+00 - r ) - r * w ) + w1

    else if ( x <= 0.18D+00 ) then

        h = x
        w1 = 0.0D+00

        r = h / ( h + two )
        t = r * r
        w = ( ( p2 * t + p1 ) * t + p0 ) / ( ( q2 * t + q1 ) * t + 1.0D+00 )
        rlog1 = two * t * ( 1.0D+00 / ( 1.0D+00 - r ) - r * w ) + w1

    else if ( x <= 0.57D+00 ) then

        h = 0.75D+00 * x - 0.25D+00
        w1 = b + h / 3.0D+00

        r = h / ( h + 2.0D+00 )
        t = r * r
        w = ( ( p2 * t + p1 ) * t + p0 ) / ( ( q2 * t + q1 ) * t + 1.0D+00 )
        rlog1 = two * t * ( 1.0D+00 / ( 1.0D+00 - r ) - r * w ) + w1

    else

        w = ( x + half ) + half
        rlog1 = x - log ( w )

    end if

    return
    end function rlog1

    subroutine student_cdf_values ( n_data, a, x, fx )

    !*****************************************************************************80
    !
    !! STUDENT_CDF_VALUES returns some values of the Student CDF.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    02 June 2001
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Milton Abramowitz, Irene Stegun,
    !    Handbook of Mathematical Functions,
    !    US Department of Commerce, 1964.
    !
    !  Parameters:
    !
    !    Input/output, integer ( kind = int32 ) N_DATA.  The user sets N_DATA to 0
    !    before the first call.  On each call, the routine increments N_DATA by 1,
    !    and returns the corresponding data; when there is no more data, the
    !    output value of N_DATA will be 0 again.
    !
    !    Output, integer ( kind = int32 ) A, real ( kind = dp ) X, the arguments of
    !    the function.
    !
    !    Output, real ( kind = dp ) FX, the value of the function.
    !
    implicit none

    integer ( kind = int32 ) :: n_data,a
    real ( kind = dp ) :: x,fx

    integer ( kind = int32 ), parameter :: n_max = 13

    integer ( kind = int32 ), save, dimension ( n_max ) :: a_vec = (/ &
        1, 2, 3, 4, &
        5, 2, 5, 2, &
        5, 2, 3, 4, &
        5 /)

    real ( kind = dp ), save, dimension ( n_max ) :: fx_vec = (/ &
        0.60D+00, 0.60D+00, 0.60D+00, 0.60D+00, &
        0.60D+00, 0.75D+00, 0.75D+00, 0.95D+00, &
        0.95D+00, 0.99D+00, 0.99D+00, 0.99D+00, &
        0.99D+00 /)
    real ( kind = dp ), save, dimension ( n_max ) :: x_vec = (/ &
        0.325D+00, 0.289D+00, 0.277D+00, 0.271D+00, &
        0.267D+00, 0.816D+00, 0.727D+00, 2.920D+00, &
        2.015D+00, 6.965D+00, 4.541D+00, 3.747D+00, &
        3.365D+00 /)

    if ( n_data < 0 ) then
        n_data = 0
    end if

    n_data = n_data + 1

    if ( n_max < n_data ) then
        n_data = 0
        a = 0
        x = 0.0D+00
        fx = 0.0D+00
    else
        a = a_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
    end if

    return
    end subroutine student_cdf_values

    function stvaln ( p )

    !*****************************************************************************80
    !
    !! STVALN provides starting values for the inverse of the normal distribution.
    !
    !  Discussion:
    !
    !    The routine returns an X for which it is approximately true that
    !      P = CUMNOR(X),
    !    that is,
    !      P = Integral ( -infinity < U <= X ) exp(-U*U/2)/sqrt(2*PI) dU.
    !
    !  Reference:
    !
    !    William Kennedy, James Gentle,
    !    Statistical Computing,
    !    Marcel Dekker, NY, 1980, page 95,
    !    QA276.4 K46
    !
    !  Parameters:
    !
    !    Input, real ( kind = dp ) P, the probability whose normal deviate
    !    is sought.
    !
    !    Output, real ( kind = dp ) STVALN, the normal deviate whose probability
    !    is approximately P.
    !
    implicit none

    real ( kind = dp ) :: p
    real ( kind = dp ) :: stvaln

    real ( kind = dp ) :: sgn,y,z

    real ( kind = dp ), parameter, dimension(0:4) :: xden = (/ &
        0.993484626060D-01, &
        0.588581570495D+00, &
        0.531103462366D+00, &
        0.103537752850D+00, &
        0.38560700634D-02 /)
    real ( kind = dp ), parameter, dimension(0:4) :: xnum = (/ &
        -0.322232431088D+00, &
        -1.000000000000D+00, &
        -0.342242088547D+00, &
        -0.204231210245D-01, &
        -0.453642210148D-04 /)

    if ( p <= 0.5D+00 ) then

        sgn = -1.0D+00
        z = p

    else

        sgn = 1.0D+00
        z = 1.0D+00 - p

    end if

    y = sqrt ( -2.0D+00 * log ( z ) )
    stvaln = y + eval_pol ( xnum, 4, y ) / eval_pol ( xden, 4, y )
    stvaln = sgn * stvaln

    return
    end function stvaln

    subroutine timestamp ( )

    !*****************************************************************************80
    !
    !! TIMESTAMP prints the current YMDHMS date as a time stamp.
    !
    !  Example:
    !
    !    31 May 2001   9:45:54.872 AM
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    06 August 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    None
    !
    implicit none

    character ( len = 8 ) :: ampm
    integer ( kind = int32 ) :: d,h,m,mm,n,s,y

    integer ( kind = int32 ) values(8)
    character ( len = 9 ), parameter, dimension(12) :: month = (/ &
        'January  ', 'February ', 'March    ', 'April    ', &
        'May      ', 'June     ', 'July     ', 'August   ', &
        'September', 'October  ', 'November ', 'December ' /)

    call date_and_time ( values = values )

    y = values(1)
    m = values(2)
    d = values(3)
    h = values(5)
    n = values(6)
    s = values(7)
    mm = values(8)

    if ( h < 12 ) then
        ampm = 'AM'
    else if ( h == 12 ) then
        if ( n == 0 .and. s == 0 ) then
            ampm = 'Noon'
        else
            ampm = 'PM'
        end if
    else
        h = h - 12
        if ( h < 12 ) then
            ampm = 'PM'
        else if ( h == 12 ) then
            if ( n == 0 .and. s == 0 ) then
                ampm = 'Midnight'
            else
                ampm = 'AM'
            end if
        end if
    end if

    write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
        d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

    return
    end subroutine timestamp

    end module cdflib

