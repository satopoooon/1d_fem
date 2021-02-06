# 対称3層スラブ構造の導波路の解析解を求める
# 参考サイト:http://www.ocw.titech.ac.jp/index.php?module=General&action=DownLoad&file=20132228841109-6-0-29.pdf&type=cal&JWC=20132228841109

const λ = 1.0
const k0 = 2pi * λ
const a = 5.0
const n1 = 1.49
const n2 = 1.48
const n = 0

module AnalyticalSolution
# TEモードの伝搬定数を求める
using NLsolve
export cal_TE_mode, nls 

function cal_TE_mode(β, k0, n1, n2, a)
    κ = sqrt((k0*n1)^2-β^2)
    γ = sqrt(β^2-(k0*n2)^2)

    return γ/κ - tan(κ*a)
end

function nls(func, params...; ini = [0.0])
    if typeof(ini) <: Number
        r = nlsolve((vout,vin)->vout[1]=func(vin[1],params...), [ini])
        v  = r.zero[1]
    else
        r = nlsolve((vout,vin)->vout .= func(vin,params...), ini)
        v = r.zero
    end
    return v, r.f_converged
end

end

using .AnalyticalSolution
using Plots

# nlsolve計算時の初期値
β_ini = 9.8* (k0*n1 - k0*n2)/10 + k0*n2 

β, jud = nls(cal_TE_mode, k0, n1, n2, a, ini = β_ini) 

if !jud
    error("β is not correct!!!")
end
print("β=$β")

κ = sqrt(k0^2*n1^2-β^2)
γ = sqrt(β^2-k0^2*n2^2)

# 係数は任意の値
Ae = 1.0

# TE 偶数モードの電解分布を計算
function Ey(x, a, κ, γ, B)

    #境界条件より決まる係数
    D = B * cos(κ*a)/exp(-γ*a)
    C = D

    if x >  a
        Ey = D * exp(-γ * x)
    elseif -a ≤ x ≤ a
        Ey = B * cos(κ * x)
    elseif x < -a
        Ey = C * exp(γ * x)
    end

    return Ey
end

# 電解分布を計算する範囲
x_list = collect(-10:0.1:10)


# 電解分布を計算
Ey_list = [Ey(x, a, κ, γ, Ae) for x in x_list]
plot(x_list, Ey_list)