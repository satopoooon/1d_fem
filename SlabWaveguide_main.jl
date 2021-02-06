module SlabWaveguide_FEM
# 1次元有限要素法でTEモードの電界分布を求める
# TEモードに関してスラブ導波路をFEMで求める
# 前提として導波路構造の媒質が等方性・非磁性(比透磁率＝1)
# 無損失（電流密度と電荷密度がJ=ρ=0）の誘電体からなるとする　　
# 光波の伝搬方向はZ,
# x方向に対して構造は不均一

# λ:波長
# k0:波数
# 2a:導波路幅
# n1:コア屈折率
# n2:クラッド屈折率
# n:自然数
const λ = 1.0
const k0 = 2pi / λ
const a = 5.0
const n1 = 1.49
const n2 = 1.48
const n = 0

export All_element, make_matrix_C

# 構造体を定義
mutable struct All_element
    N::Int64 # nodeの数
    seg::Array{Array{Float64,1},1} # 線要素を格納する配列
end

# 多重ディスパッチにより,引数の型が「Array{Float64,1}」のときは下記関数が実行される
function All_element(x::Array{Float64,1})
    N = length(x)
    seg = make_segments(x)
    
    # 多重ディスパッチにより構造体のAll_elementが実行される
    return All_element(N, seg)

end

function make_segments(vec_x)
    N = length(vec_x)                                           # vec_xの要素数
    vec_segments = Array{Array{Float64,1},1}(undef,N-1)         # r_i,r_i+1 を配列に格納する．配列数はelementの数(vec_x-1になる)分 
    for i=1:N-1                                                 
        vec_segments[i] = [vec_x[i],vec_x[i+1]]                 
    end
    return vec_segments
end

# All_element構造体のメソッド
function make_matrix_C(self::All_element)
    seg = self.seg
    N = self.N
    C = zeros(Float64,N, N)

    for e=1:N-1

        tmp_element = Element(seg[e])

        Ce = tmp_element.Ce

        C[e, e]     = C[e, e] + Ce[1,1]
        C[e+1, e]   = C[e+1, e] + Ce[2,1]
        C[e, e+1]   = C[e, e+1] + Ce[1,2]
        C[e+1, e+1] = C[e+1, e+1] + Ce[2,2]

    end

    return C

end

struct Element
    seg::Array{Float64,1}
    len::Float64
    Ce::Array{Float64,2}

    # コンストラクタ(引数の数が構造体と同じ場合は下記のように書く)
    # 尚、All_Element構造体のコンストラクタは、多重ディスパッチにより配列を引数とするAll_Element関数を実行して作成したので作り方が異なる．
    function Element(seg)
        len = calc_length(seg)
        Ce = C_element(len, seg)
        
        return new(seg, len, Ce)
    end

    function calc_length(vec_segments)
        edge_n = length(vec_segments) - 1 # エッジ数=ノード数 - 1
        len = zeros(Float64,edge_n)

        len = vec_segments[2]-vec_segments[1]

        return len
    end

    # 各nodeでのCeを計算する．
    function C_element(len, vec_segments)
        
        # 屈折率を座標により変える 
        if -a ≤ vec_segments[1] &&  vec_segments[2] ≤　a  
            n_e = 1.49
        else
            n_e = 1.48
        end

        len = len
        #n_e = 0
        #n_e = 999999999999
        c11 = -6/len^2 + n_e^2 * k0^2
        #println(c11)
        c22 = c11
        c12 = 6/len^2
        c21 = c12        
        
        return [c11 c12; c21 c22]
    end
    
end

end

module AnalyticalSolution
# 対称3層スラブ構造の導波路の解析解を求める
# 参考サイト:http://www.ocw.titech.ac.jp/index.php?module=General&action=DownLoad&file=20132228841109-6-0-29.pdf&type=cal&JWC=20132228841109
# TEモードの伝搬定数を求める
using NLsolve
export cal_TE_mode, nls, Ey

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

# TE 偶数モードの電解分布を計算
function Ey(x, a, κ, γ)
    
    B =1
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

end

using .SlabWaveguide_FEM
using .AnalyticalSolution
using Plots
using LinearAlgebra

function main_AnalySol()

    λ = 1.0
    k0 = 2pi * λ
    a = 5.0
    n1 = 1.49
    n2 = 1.48
    
    #using .AnalyticalSolution
    #using LinearAlgebra
    #using Plots
    
    # nlsolve計算時の初期値
    β_ini = 9.8* (k0*n1 - k0*n2)/10 + k0*n2 

    β, jud = nls(cal_TE_mode, k0, n1, n2, a, ini = β_ini) 

    if !jud
        error("β is not correct!!!")
    end
    print("β=$β")

    κ = sqrt(k0^2*n1^2-β^2)
    γ = sqrt(β^2-k0^2*n2^2)

    # 電解分布を計算する範囲
    x_list = collect(-20:0.1:20)

    # 電解分布を計算
    Ey_list = [Ey(x, a, κ, γ) for x in x_list]
    #println(Ey_list)

    # FEMと比較するために電解ベクトルのnormが１になるように正規化
    Ey_list = Ey_list/norm(Ey_list)

    #plot(x_list[2:end-1], Ey_list[2:end-1])
    return x_list[2:end-1], Ey_list[2:end-1]
end

function main_FEM()
# 有限要素法を実施する座標を定義
    x_min = -20
    x_max = 20
    step = 0.1                                                 
    x = collect(x_min:step:x_max)

    ele = All_element(x)
    C = make_matrix_C(ele)

    C = C[2:end-1,2:end-1]

    # 固有値及び固有ベクトルを求める
    β, vec = eigen(C)
    x = x[2:end-1]
    return x, vec[:,end]
end

# 電解成分を計算
x1, AnalySol = main_AnalySol()
x2, FEMSol = main_FEM()

plot(x2, FEMSol)
plot!(x1, AnalySol)