import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
 
# mk system
def func(state,t, m, k):
    x1, x2 = state
    dx2dt = - (k/m) * x1
    return [x2, dx2dt]
 
m = 0.5                         # 質量[kg]
k = 1000                        # 剛性[N/m]
state0 = [0.0, 0.1]             # 初期値[x0, v0]
 
t0 = 0                          # 初期時間[s]
tf = 1                          # 終了時間[s]
dt = 0.005                      # 時間刻み[s]
t = np.arange(t0, tf+dt, dt)    # 時間軸配列
 
# odeintを使った解法
# x1, x2の解を求めているので、sol[0]が変位x1, sol[1]が速度x2となる
sol = odeint(func, state0, t, args=(m,k))
 
 
omega = np.sqrt(k/m)
theory = (state0[0] * omega * np.cos(omega * t)) +\
         ((state0[1]/omega) * np.sin(omega * t))
 
# ここからグラフ描画
# フォントの種類とサイズを設定する。
plt.rcParams['font.size'] = 14
plt.rcParams['font.family'] = 'Times New Roman'
 
# 目盛を内側にする。
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
 
# グラフの上下左右に目盛線を付ける。
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
 
# 軸のラベルを設定する。
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Displacement [m]')
 
# データの範囲と刻み目盛を明示する。
ax1.set_xticks(np.arange(0, 2, 0.2))
ax1.set_yticks(np.arange(-1, 1, 0.001))
ax1.set_xlim(0, 1)
ax1.set_ylim(-0.005, 0.005)
 
# データプロット
ax1.plot(t, sol[:,0], label='Python odeint [x0=0, v0=0.1]', c='b', marker='o', linestyle='None')
ax1.plot(t, theory, label='Theory', c='r')
 
fig.tight_layout()
plt.legend(loc='upper left')
 
# グラフを表示する。
plt.show()
plt.close()