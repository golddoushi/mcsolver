# mcsolver
A user friendly tools implementing Monte Carlo simulations to estimate Curie/Neel temperature

Original version contributor： Dr. Liang Liu* 1.Shenzheng University 2.Shandong University
Email: liangliu@mail.sdu.edu.cn

You can download the packed .exe (only tested in Windows 10 platform) from the following link. Wish it can find something helpful for you. And if it was used for publication, please cite:
[1] Magnetic switches via electric field in BN nanoribbons. Applied Surface Science 480(2019)

网盘下载地址
链接：https://pan.baidu.com/s/1EaDqOOdB7AP9WXrwEIEaxQ
提取码：52ze

安装方法：
无需安装

使用方法：
打开软件（打开较慢约10s），从上至下依次填写参数，然后点击startMC即可。如有帮助请引用论文，谢谢。

1.填写/修改晶格基矢

2.增查改删轨道信息，注意坐标为分数坐标。Ani项后面是三个方向的single-ion anisotropy（Ising模型无此项，xy模型需要设定前两个也就是Dz，Dx），注意此处以及下面所有的能量的单位都是K，与meV的换算见下。

3.增查改删交换作用（bond），包括xyz三个方向的交换强度（如果用Ising模型则只需要设定jz，xy模型需要设定前两个）、交换链接的两个轨道的id（在上面一步中定义了）、交换跨越的晶格矢量。
点选一个列表中的交换，结构预览中的键会变成粗黄线，看看是否与预料中的一样。
轨道与交换作用修改之后，或者点选交换作用的列表，结构预览就会更新，可以用鼠标左键拖动可翻转，右键拖动可放缩，多角度检查交换构型。

4.设定其他参数，包括温度始末点以及总的温度插值点、nthermal：热化（达到热平衡）所需的MC步、nsweep：热化后的统计次数、模型、算法（暂只支持metroplis局域更新，与wolff区块更新，即将加入Swensden-Wang算法、continuous time 算法与order conserving算法）

5.设定并行线程数

6.(可选)点击save按钮保存当前设置

7.点击startMC按钮

8.右侧出图之后就计算完成了。在软件所在目录有一个result.txt记录了平均自旋、磁化率、能量、比热、Binder cumulate U4等信息。如果并行计算温度的次序可能是错乱的，但是每行的对应是正确的。

9.(可选)可以点击load按钮载入已有设置。sample文件夹中为CrI3的Heisenberg模型设置，可参考。设置文件为txt格式，可直接用记事本进行修改。