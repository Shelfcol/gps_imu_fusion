### Forked from https://github.com/zm0612/eskf-gps-imu-fusion

# ESKF融合IMU与GPS数据
结果如下：

![融合IMU数据之后的GPS轨迹效果](/data/raw_data1/results/trajectory.png)
![融合IMU数据之后的GPS轨迹效果](/data/raw_data1/results/xyz_view.png)

绿色轨迹：ground truth
蓝色轨迹：fuse IMU and GPS
红色轨迹：GPS

本人的实现方法请参考博客[《Quaternion kinematics for error state kalman filter实现GPS+IMU融合，（附源码）》](https://blog.csdn.net/qq_38650944/article/details/123580686)

欢迎大家交流呀！
## 1.  依赖库

Eigen

```shell
sudo apt-get install libeigen3-dev
```

Yaml

```shell
sudo apt-get install libyaml-cpp-dev
```

## 2. 编译

```shell
cd eskf-gps-imu-fusion
mkdir build
cd build
cmake ..
make 
```

## 3. 运行

```shell
cd eskf-gps-imu-fusion/build
./gps_imu_fusion
```

## 4.轨迹显示

执行完`./gps_imu_fusion`会生成轨迹文件

```shell
cd eskf-gps-imu-fusion/data
evo_traj kitti fused.txt gt.txt measured.txt -p
```

> 需要安装evo，可以参考博客中的介绍：https://blog.csdn.net/u011341856/article/details/104594392?spm=1001.2014.3001.5501


