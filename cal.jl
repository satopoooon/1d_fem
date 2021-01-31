# 対称3層スラブ構造の導波路の解析解を求める
# 参考サイト:http://www.ocw.titech.ac.jp/index.php?module=General&action=DownLoad&file=20132228841109-6-0-29.pdf&type=cal&JWC=20132228841109


# TEモードのEven modeの伝搬定数を求める特性方程式
using NLsolve
using Plots

#function cal_TE_mode(β, n1, n2, k0, a, n)
# function cal_TE_mode(b, V)
#     # β:伝搬定数
#     # n1:コア屈折率
#     # n2:クラッド屈折率
#     # k0:波数
#     # 2a:コア幅
#     # n:自然数

#     # # b:規格化伝搬定数
#     #b = ( (β/k0)^2 - n2^2 ) / ( n1^2 - n2^2 ) 
#     # # Δ:比屈折率差Δ
#     # Δ = (n1^2 - n2^2)/2n1^2
#     # # V:Vパラメータ
#     # V = k0 * n1 * a * sqrt(2Δ)

#     return 1/real(sqrt(complex(1-b))) * (atan(real(sqrt(complex(b/(1-b)))) + (n*pi)/2)) - V

# end

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

λ = 1.0
k0 = 2pi * λ
a = 5.0
n1 = 1.49
n2 = 1.48
n = 0

# b:規格化伝搬定数
#b = ( (β/k0)^2 - n2^2 ) / ( n1^2 - n2^2 ) 
# Δ:比屈折率差Δ
# Δ = (n1^2 - n2^2)/2n1^2
# # V:Vパラメータ
# V_ini = k0 * n1 * a * sqrt(2Δ)
# println(V_ini)

#b_cal, jud = nls(cal_TE_mode, V_ini, ini = 0.9) 
# println(b_cal)
# if jud == false
#     println("cal_b is not correct!!!") 
# else
#     β = sqrt(b_cal*(n1^2 - n2^2)*k0^2 + n2^2)
# end

#κ

print("cal start")
β_ini = 9.8* (k0*n1 - k0*n2)/10 + k0*n2 
β, jud = nls(cal_TE_mode, k0, n1, n2, a, ini = β_ini) 
print(β)
print(jud)

κ = sqrt(k0^2*n1^2-β^2)
γ = sqrt(β^2-k0^2*n2^2)

Ae = 1.0
# TE 偶数モードの電解
# コア部
function Ey(x, a, κ, γ, B)

    D = B * cos(κ*a)/exp(-γ*a)
    C = D
    if x >  a
        Ey = D * exp(-γ * x)
        #Ey = 0
    elseif -a ≤ x ≤ a
        Ey = B * cos(κ * x)
    elseif x < -a
        Ey = C * exp(γ * x)
        #Ey = 0 
    end

    return Ey
end

x_list = collect(-10:0.1:10)
#print(x_list)

Ey_list = [Ey(x, a, κ, γ, Ae) for x in x_list]
plot(x_list, Ey_list)

# v = cal_TE_mode(0.1)
# println(v)
# v = cal_TE_mode(0.5)
# println(v)
# v = cal_TE_mode(0.9)
# println(v)
# v = cal_TE_mode(0.95)
# println(v)


#ret = cal_TE_mode(0.8, n1, n2, k0, a, n)
