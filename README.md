# mcsolver
A user friendly tools using Monte Carlo simulations for estimation of Curie temperature

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
打开软件（打开较慢约10s），从上至下依次填写参数，然后点击submit即可。如有帮助请引用论文，谢谢。

1.填写/修改晶格基矢

2.增查改删轨道信息，注意坐标为分数坐标

3.增查改删交换作用（Bond），包括xyz三个方向的交换强度（目前只需要设定Jz）、交换链接的两个轨道的ID、交换跨越的晶格矢量。
   点选一个列表中的交换，结构预览中的键会变成粗黄线，看看是否与符合设想，不要设置冗余的交换（即不能有线条首尾完全重合）。
   修改轨道、交换作用，或者点选交换作用的列表，结构预览都会更新，可以用鼠标翻转结构，多角度检查轨道、耦合是否正确。

4.设定其他参数，包括温度始末点以及温度取样数目；nthermal：热化所需次数；nsweep：热化后的统计次数；模型（暂时只支持Ising，即将
   增加XY模型，Heisenberg模型等。不久之后还会加入部分量子蒙特卡洛模型）；算法暂只支持Metroplis局域更新（适合模拟远离相变的温
   区）与Wolff区块更新（适合相变点附近的温区），以后会加入Sweden-Wang、continuous time与Order conserving等算法。

5.设定并行线程数cores

6.点击submit按钮

7.右侧出图（平均净自旋对各温度的散点图）之后就计算完成了。在软件所在目录有一个result.txt记录了平均自旋、磁化率、能量、比热等信息。
   如果采用并行计算（即ncores>1）温度的次序可能是错乱的，但是每行的对应是正确的。
