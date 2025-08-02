using Base.GMP.MPZ: sizeinbase, mpz_t, bitcnt_t
using Base.GMP: CulongMax, ClongMax, libgmp

if Clong == Int32
    const LargeInt = Union{Int64, Int128}
    const LargeUInt = Union{UInt64, UInt128}
else
    const LargeInt = Int128
    const LargeUInt = UInt128
end

gmpz(op::Symbol) = (Symbol(:__g, op), libgmp)

# mpz_t -> BigInt
# bitcnt_t -> bitcnt_t
# Culong -> CulongMax
# Clong -> ClongMax
function cnv(op::Symbol)::Symbol
    if op === :mpz_t
        return :BigInt
    elseif op === :bitcnt_t
        return :bitcnt_t
    else
        Symbol(op, :Max)
    end
end


# addmul, submul
for (fname, gmpname, ytype) in [
    (:addmul_mpz!, :mpz_addmul,    :mpz_t ),
    (:addmul_ul!, :mpz_addmul_ui, :Culong),
    (:submul_mpz!, :mpz_submul,    :mpz_t ),
    (:submul_ul!, :mpz_submul_ui, :Culong),]

    @eval begin
        function $fname(z::BigInt, x::BigInt, y::$(cnv(ytype)))
            ccall($(gmpz(gmpname)), Cvoid, (mpz_t, mpz_t, $ytype), z, x, y)
        end
    end
end


# isdivisible
for (fname, gmpname, ytype) in [
    (:isdivisible_mpz,  :mpz_divisible_p,      :mpz_t   ),
    (:isdivisible_ul,  :mpz_divisible_ui_p,   :Culong  ),
    (:isdivisible_2exp, :mpz_divisible_2exp_p, :bitcnt_t),]

    @eval begin
        function $fname(n::BigInt, d::$(cnv(ytype)))::Bool
            !iszero(ccall($(gmpz(gmpname)), Cint, (mpz_t, $ytype), n, d))
        end
    end
end


# divexact
for (fname, gmpname, ytype) in [
    (:divexact_mpz, :mpz_divexact,    :mpz_t ),
    (:divexact_ul, :mpz_divexact_ui, :Culong),]

    @eval begin
        function $fname(n::BigInt, d::$(cnv(ytype)))::BigInt
            z = BigInt()
            ccall($(gmpz(gmpname)), Cvoid, (mpz_t, mpz_t, $ytype), z, n, d)
            return z
        end
    end
end


# mod_ul
mod_ul(a::BigInt, b::CulongMax)::Culong = ccall((:__gmpz_fdiv_ui, libgmp), Culong, (mpz_t, Culong), a, b)


# iscongruent
for (fname, gmpname, ytype1, ytype2) in [
    (:iscongruent_mpz,  :mpz_congruent_p,      :mpz_t,  :mpz_t   ),
    (:iscongruent_ul,  :mpz_congruent_ui_p,   :Culong, :Culong  ),
    (:iscongruent_2exp, :mpz_congruent_2exp_p, :mpz_t,  :bitcnt_t),]

    @eval begin
        function $fname(n::BigInt,c::$(cnv(ytype1)), d::$(cnv(ytype2)))::Bool
            !iszero(ccall($(gmpz(gmpname)), Cint, (mpz_t, $ytype1, $ytype2), n, c, d))
        end
    end
end


# powm_ui
function powm_ul(base::BigInt,exp::CulongMax,mod::BigInt)
    rop = BigInt()
    ccall(gmpz(:mpz_powm_ui), Cvoid, (mpz_t, mpz_t, Culong, mpz_t), rop, base, exp, mod)
    return rop
end


# ui_pow_ui
# If you are computing a power x^y between integers x,y that fits into Culong,
# you do not need to set the result as BigInt(x)^y, even if the result may exceed Int128.
# just pow_ul_ul(x,y).
# もしあなたがCulongに収まる程度の整数x,y同士のべき乗x^yを計算する場合、その結果がInt128を超える可能性があってもBigInt(x)^yとする必要はなく、pow_ul_ul(x,y)とすればいいです。
"""
    pow_ul_ul(base, exp)

Return `base^exp` as a `BigInt`. Wrapper for `__gmpz_ui_pow_ui`.

Useful for fast small exponentiation without promotion to `BigInt` until the final step.
"""
function pow_ul_ul(base::CulongMax, exp::CulongMax)::BigInt
    rop = BigInt()
    ccall(gmpz(:mpz_ui_pow_ui), Cvoid, (mpz_t, Culong, Culong), rop, base, exp)
    return rop
end


# root
function iroot_ul(x::BigInt, n::CulongMax)::BigInt
    z = BigInt()
    ccall(gmpz(:mpz_root), Cint, (mpz_t, mpz_t, Culong), z, x, n)
    return z
end
function rootrem_ul(x::BigInt, n::CulongMax)::BigInt
    rem = BigInt()
    root = BigInt()
    ccall(gmpz(:mpz_rootrem), Cvoid, (mpz_t, mpz_t, mpz_t, Culong), root, rem, x, n)
    return (root, rem)
end
function sqrtrem(x::BigInt)::BigInt
    rem = BigInt()
    root = BigInt()
    ccall(gmpz(:mpz_sqrtrem), Cvoid, (mpz_t, mpz_t, mpz_t), root, rem, x)
    return (root, rem)
end


# perfect power and square
for (fname, gmpname) in [
    (:isperfectpower,  :mpz_perfect_power_p) ,
    (:isperfectsquare, :mpz_perfect_square_p),]

    @eval begin
        function $fname(n::BigInt)::Bool
            !iszero(ccall($(gmpz(gmpname)), Cint, (mpz_t,), n))
        end
    end
end


# PRIME-related is already implemented in Primes.jl


# gcd and lcm
function gcd_ul(op1::BigInt, op2::CulongMax)::Culong
    iszero(op2) && throw(DomainError(op2, "0 cannot be specified for op2 in gcd_ul()"))
    return ccall(gmpz(:mpz_gcd_ui), Culong, (Ptr{Cvoid}, mpz_t, Culong), C_NULL, op1, op2)
end

function lcm_ul(op1::BigInt, op2::CulongMax)::BigInt
    rop = BigInt()
    ccall(gmpz(:mpz_lcm_ui), Cvoid, (mpz_t, mpz_t, Culong), rop, op1, op2)
    return rop
end


# jacobi and kronecker
for (fname, gmpname, ytype1, ytype2) in [
    (:jacobi,   :mpz_jacobi,   :mpz_t, :mpz_t),
    (:legendre, :mpz_legendre, :mpz_t, :mpz_t),
    (:kronecker,    :mpz_kronecker,    :mpz_t,  :mpz_t ),
    (:kronecker_sl, :mpz_kronecker_si, :mpz_t,  :Clong ),
    (:kronecker_ul, :mpz_kronecker_ui, :mpz_t,  :Culong),
    (:kronecker_sl, :mpz_si_kronecker, :Clong,  :mpz_t ),
    (:kronecker_ul, :mpz_ui_kronecker, :Culong, :mpz_t ),]

    @eval begin
        function $fname(n::$(cnv(ytype1)),k::$(cnv(ytype2)))::Cint
            return ccall($(gmpz(gmpname)), Cint, ($ytype1, $ytype2), n, k)
        end
    end
end