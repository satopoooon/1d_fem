# TEモードに関してスラブ導波路をFEMで求める
# 前提として導波路構造の媒質が等方性・非磁性(比透磁率＝1)
# 無損失（電流密度と電荷密度がJ=ρ=0）の誘電体からなるとする　　
# 光波の伝搬方向はZ,
# x方向に対して構造は不均一

# 1次元有限要素法でTEモードの電界分布を求める
module SlabWaveguide_FEM

# λ:波長
# k0:波数
# 2a:導波路幅
# n1:コア屈折率
# n2:クラッド屈折率
# n:自然数
const λ = 1*10^-6
const k0 = 2pi / λ
const a = 6.0
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
        #n_e = 0
        #n_e = 999999999999
        c11 = -6/len^2 + n_e^2 * k0^2
        #println(c11)
        c22 = c11
        #println(c22)
        c12 = 6/len^2
        #println(C12)
        c21 = c12        
        
        return [c11 c12; c21 c22]
    end
    
end

end

using Plots
using .SlabWaveguide_FEM
using LinearAlgebra

# 有限要素法を実施する座標を定義
x_min = -10
x_max = 10
step = 0.5                                                 
x = collect(x_min:step:x_max)

ele = All_element(x)
C = make_matrix_C(ele)
#println(C)
# const λ = 1.0
# const k0 = 2pi * λ

# len = 0.1
# n_e = 1.49
# c11 = -6/len^2 + n_e^2 * k0^2
# c22 = c11
# c12 = 6/len^2
# c21 = c12    
   
# C = [c11 c12; c21 c22]

# 固有値及び固有ベクトルを求める
β, vec = eigen(C)
#print(β)
#print(vec)
# 固有値最大のときの固有ベクトルが基本モードでの電解分布となる
#plt = plot(x, real(vec[:,end]))

#println(vec[:,end])
#scatter(x, vec[:,end])


vec1 = vec[:,end] 
scatter(x, vec1)
#print(vec1)

# vec2 = [(if x < 1.0e-10 0 else x end) for x in vec1]
# println(typeof(vec2[1]))

# open("C:/Users/satoshi/julia/1d_femtest.txt", "a") do f
#     for i in vec1
#         write(f, i)
#     end
# end
