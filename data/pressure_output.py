# ライブラリのインポート
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as cm

# "flow3.csv"(pressureデータの入ったファイル)を読み込む
flow3 = pd.read_csv('flow3.csv', header=None)

# x, y方向それぞれの格子番号
x = range(0, 401)
y = range(0, 201)

# 等高線を描画するためにmeshgridに変換
X, Y = np.meshgrid(x, y)

# flow_lsはpressureデータのみを抽出したもの(indexを取り除いた)
flow_ls = list(flow3[0])

# アニメーションの出力
## imgsに各々の出力ステップでのimgを入れていく
fig = plt.figure(figsize=(10, 6))
imgs = []
for n in range(20):
    p = []
    flow = flow_ls[401 * 201 * n: 401 * 201 * (n + 1)]
    for i in range(401):
        tmp = flow[201 * i: 201 * (i + 1)]
        p.append(tmp)
    p = np.array(p).T
    img = plt.contour(X, Y, p, 100, cmap=cm.bwr)
    img = img.collections
    imgs.append(img)

plt.title('Pressure distribution')
plt.xlabel('grid number in x axis')
plt.ylabel('grid number in y axis')
pp=plt.colorbar (orientation="vertical")

anim = animation.ArtistAnimation(fig, imgs, interval=160)
anim.save('pressure.gif', 'imagemagick')

plt.show()
