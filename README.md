# 2D RVE Model Generator for ABAQUS

[English](#english) | [中文](#chinese)

<a name="english"></a>

## 📋 Overview

A powerful Python script for automated generation of 2D Representative Volume Element (RVE) models in ABAQUS with controllable fiber volume fraction. This tool is designed for composite material simulation and micromechanical analysis.

### ✨ Key Features

- **Controllable Fiber Distribution**: Automatic generation of 2D RVE geometric models based on user-defined fiber volume fraction (Vf)
- **Advanced Placement Algorithm**: Hybrid algorithm combining RSA (Random Sequential Adsorption) seeding, anchored relaxation, and forced correction for optimal fiber distribution
- **Minimum Distance Guarantee**: Ensures minimum spacing between fibers with error termination if requirements cannot be met
- **Complete Material Definition**: Automatic creation of materials (matrix + fiber + interface) with section assignment and mesh generation
- **Periodic Boundary Conditions (PBCs)**: Implements XY direction periodicity for accurate RVE behavior
- **Cohesive Elements**: Automatically inserts 2D cohesive elements (COH2D4) at fiber-matrix interfaces
- **Accurate Phase Recognition**: Geometric location algorithm for precise identification of fibers and matrix (including corner fragments)
- **Data Export**: Exports fiber center coordinates to CSV files for verification and analysis
- **Distance Verification**: Validates inter-fiber spacing with detailed statistical output

## 🎯 Applications

- Composite material micromechanics simulation
- Fiber-reinforced composite analysis
- Interfacial debonding studies
- Homogenization analysis
- Material property prediction
- Damage and failure mechanism investigation

## 🔧 Requirements

### Software Requirements
- **ABAQUS**: Version 2023 or compatible
- **Python**: 2.7 (ABAQUS built-in)

### Python Modules
All required modules are included in ABAQUS installation:
- `abaqus`, `abaqusConstants`, `caeModules`
- Standard library: `math`, `random`, `time`, `os`

## 📥 Installation

1. Clone this repository:
```bash
git clone https://github.com/ZPL-03/2D_RVE_Model.git
cd 2D_RVE_Model
```

2. Ensure ABAQUS is properly installed and accessible from command line

## 🚀 Quick Start

### Basic Usage

1. **Edit Parameters** (Lines 1301-1494 in `2D_RVE_Model.py`):
```python
# Basic geometry parameters
RVE_SIZE = [1.0, 1.0]          # RVE dimensions [width, height] in mm
FIBER_RADIUS = 0.05             # Fiber radius in mm
TARGET_VF = 0.50                # Target fiber volume fraction (0.0-1.0)
MIN_DIST_FACTOR = 0.01          # Minimum distance factor
```

2. **Run in ABAQUS CAE**:
   - Open ABAQUS CAE
   - File → Run Script → Select `2D_RVE_Model.py`
   
   Or from command line:
   ```bash
   abaqus cae noGUI=2D_RVE_Model.py
   ```

### Advanced Configuration

#### Distribution Control
```python
# RSA seeding ratio - controls fiber distribution pattern
RSA_SEEDING_RATIO = 0.9  # 0.0-1.0
# High values (0.8-1.0): Fast, uniform distribution
# Low values (0.0-0.3): Realistic clustering, slower generation
# Medium (0.4-0.7): Balanced approach
```

#### Mesh Settings
```python
GLOBAL_SEED_SIZE = 0.01  # Global mesh seed size in mm
```

#### Material Parameters
```python
# Fiber properties (MPa)
FIBER_E = 15000.0          # Elastic modulus
FIBER_NU = 0.3             # Poisson's ratio

# Matrix properties (MPa)
MATRIX_E = 3170.0          # Elastic modulus
MATRIX_NU = 0.35           # Poisson's ratio

# Cohesive interface properties
COHESIVE_K_NN = 1e8        # Normal stiffness (N/mm³)
COHESIVE_T_N = 44.0        # Normal strength (MPa)
COHESIVE_GIC = 0.001       # Mode I fracture energy (N/mm)
```

## 📊 Output Files

### ABAQUS Model
- Model name: `2D-RVE-Vf-XX` (XX = volume fraction percentage)
- Includes:
  - Complete geometry with fibers and matrix
  - Material definitions
  - Cohesive interface elements
  - Periodic boundary conditions
  - Meshed assembly ready for analysis

### CSV Export (Optional)
- **Filename**: `FiberCenters_VfXX_YYYYMMDD_HHMMSS.csv`
- **Content**: 
  - Fiber center coordinates
  - RVE dimensions
  - Volume fraction information
  - Timestamp and metadata

**Sample CSV Format**:
```csv
# RVE Fiber Center Coordinates
# Generated: 2025-10-20 14:30:00
# RVE Size (Width x Height): 1.000000 x 1.000000
# Fiber Radius: 0.050000
# Target Volume Fraction: 0.5000
# Total Fiber Count (effective in RVE): 127
Fiber_ID,X_Coordinate,Y_Coordinate
1,0.12345678,0.23456789
2,0.34567890,0.45678901
...
```

## 📖 Algorithm Details

### Fiber Placement Strategy

The script employs a sophisticated three-phase algorithm:

1. **RSA Seeding Phase**: 
   - Generates initial fiber positions using Random Sequential Adsorption
   - Controlled by `RSA_SEEDING_RATIO` parameter
   - Ensures no overlap during initial placement

2. **Anchored Relaxation Phase**:
   - Refines fiber positions to achieve target volume fraction
   - Maintains minimum distance constraints
   - Iteratively adjusts positions for optimal packing

3. **Forced Correction Phase**:
   - Final adjustment to meet exact volume fraction requirements
   - Handles boundary conditions and periodic constraints

### Periodic Boundary Conditions

The implementation includes:
- Automatic node pairing on opposite boundaries
- Equation constraints for displacement continuity
- Reference points for controlling periodic deformation
- Support for both X and Y directions

### Cohesive Interface Implementation

- Automatic detection of fiber-matrix interfaces
- Insertion of COH2D4 elements
- Support for mixed-mode fracture criteria (BK criterion)
- Viscous stabilization for numerical stability

## 🎨 Customization Guide

### Adjusting Volume Fraction
```python
TARGET_VF = 0.30  # For 30% fiber volume fraction
```
**Note**: Maximum achievable Vf depends on fiber radius and RVE size. Script will report errors if target cannot be met.

### Changing Distribution Pattern
```python
# For uniform, production-ready models
RSA_SEEDING_RATIO = 0.9

# For realistic, clustered distributions
RSA_SEEDING_RATIO = 0.1

# For balanced approach
RSA_SEEDING_RATIO = 0.5
```

### Modifying Mesh Density
```python
# Finer mesh (more accurate, slower)
GLOBAL_SEED_SIZE = 0.005

# Coarser mesh (faster, less accurate)
GLOBAL_SEED_SIZE = 0.02
```

## ⚠️ Important Notes

### Limitations
- Maximum fiber volume fraction depends on geometric constraints
- Cohesive element insertion may increase computational cost
- Periodic boundary conditions require careful load application

### Best Practices
1. Start with lower volume fractions (Vf < 0.6) for initial testing
2. Adjust `MIN_DIST_FACTOR` based on mesh density requirements
3. Verify fiber spacing using the exported CSV coordinates
4. Check mesh quality before running analyses
5. Use appropriate material parameters for your specific composite system

### Troubleshooting

**Problem**: Script reports "Cannot achieve target volume fraction"
- **Solution**: Decrease `TARGET_VF` or increase `RVE_SIZE` or decrease `FIBER_RADIUS`

**Problem**: Mesh generation fails
- **Solution**: Increase `GLOBAL_SEED_SIZE` or adjust geometry

**Problem**: Analysis convergence issues
- **Solution**: Check cohesive parameters, adjust `COHESIVE_STAB_COEFF`

## 📝 Citation

If you use this code in your research, please cite:

```bibtex
@software{2d_rve_model_2025,
  author = {Liu, Zhengpeng},
  title = {2D RVE Model Generator for ABAQUS},
  year = {2025},
  url = {https://github.com/ZPL-03/2D_RVE_Model},
  version = {1.0}
}
```

## 👥 Author

**Liu Zhengpeng (刘正鹏)**
- GitHub: [@小盆i](https://github.com/ZPL-03)
- Email: 1370872708@qq.com / Zhengpeng0105@gmail.com
- Technical Blog: CSDN/知乎 @小盆i

## 🤝 Contributing

Contributions are welcome! Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct and the process for submitting pull requests.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- Inspired by composite micromechanics research community
- Thanks to all contributors and users for feedback and improvements

## 📮 Contact & Support

- **Issues**: Please report bugs and request features via [GitHub Issues](https://github.com/ZPL-03/2D_RVE_Model/issues)
- **Discussions**: Join discussions in [GitHub Discussions](https://github.com/ZPL-03/2D_RVE_Model/discussions)
- **Email**: For private inquiries, contact Zhengpeng0105@gmail.com

---

<a name="chinese"></a>

# ABAQUS二维RVE建模工具

## 📋 项目概述

这是一个功能强大的Python脚本，用于在ABAQUS中自动生成具有可控纤维体积分数的二维代表性体积单元（RVE）模型。该工具专为复合材料仿真和微观力学分析设计。

### ✨ 核心功能

- **可控纤维分布**：根据用户设定的纤维体积分数(Vf)自动生成二维RVE几何模型
- **先进的布置算法**：采用混合算法（RSA播种 + 锚定松弛 + 强制校正）实现最优纤维分布
- **最小距离保证**：确保纤维间的最小间距，无法满足时报错终止
- **完整材料定义**：自动创建材料（基体+纤维+界面）、赋予截面、划分网格
- **周期性边界条件**：实现XY方向周期性边界条件，准确模拟RVE行为
- **内聚力单元**：在纤维-基体界面自动插入二维内聚力单元(COH2D4)
- **精确相识别**：利用几何位置算法准确识别纤维和基体（包括边角碎片）
- **数据导出**：导出纤维中心坐标为CSV文件，便于验证和后续分析
- **距离验证**：验证纤维间距并输出详细统计信息

## 🎯 应用场景

- 复合材料微观力学仿真
- 纤维增强复合材料分析
- 界面脱粘研究
- 均匀化分析
- 材料性能预测
- 损伤失效机理研究

## 🔧 系统要求

### 软件要求
- **ABAQUS**：2023版本或兼容版本
- **Python**：2.7（ABAQUS内置）

### Python模块
所有必需模块均包含在ABAQUS安装中：
- `abaqus`, `abaqusConstants`, `caeModules`
- 标准库：`math`, `random`, `time`, `os`

## 📥 安装说明

1. 克隆仓库：
```bash
git clone https://github.com/ZPL-03/2D_RVE_Model.git
cd 2D_RVE_Model
```

2. 确保ABAQUS已正确安装并可从命令行访问

## 🚀 快速开始

### 基本使用

1. **编辑参数**（在`2D_RVE_Model.py`第1301-1494行）：
```python
# 基本几何参数
RVE_SIZE = [1.0, 1.0]          # RVE尺寸 [宽度, 高度]，单位mm
FIBER_RADIUS = 0.05             # 纤维半径，单位mm
TARGET_VF = 0.50                # 目标纤维体积分数 (0.0-1.0)
MIN_DIST_FACTOR = 0.01          # 最小距离系数
```

2. **在ABAQUS CAE中运行**：
   - 打开ABAQUS CAE
   - 文件 → 运行脚本 → 选择 `2D_RVE_Model.py`
   
   或从命令行运行：
   ```bash
   abaqus cae noGUI=2D_RVE_Model.py
   ```

### 高级配置

#### 分布控制
```python
# RSA播种比例 - 控制纤维分布模式
RSA_SEEDING_RATIO = 0.9  # 0.0-1.0
# 高值 (0.8-1.0)：快速、均匀分布
# 低值 (0.0-0.3)：真实团簇效应，生成较慢
# 中等 (0.4-0.7)：平衡方案
```

#### 网格设置
```python
GLOBAL_SEED_SIZE = 0.01  # 全局网格种子尺寸，单位mm
```

#### 材料参数
```python
# 纤维性能 (MPa)
FIBER_E = 15000.0          # 弹性模量
FIBER_NU = 0.3             # 泊松比

# 基体性能 (MPa)
MATRIX_E = 3170.0          # 弹性模量
MATRIX_NU = 0.35           # 泊松比

# 内聚力界面性能
COHESIVE_K_NN = 1e8        # 法向刚度 (N/mm³)
COHESIVE_T_N = 44.0        # 法向强度 (MPa)
COHESIVE_GIC = 0.001       # I型断裂能 (N/mm)
```

## 📊 输出文件

### ABAQUS模型
- 模型名称：`2D-RVE-Vf-XX`（XX为体积分数百分比）
- 包含内容：
  - 完整的纤维和基体几何
  - 材料定义
  - 内聚力界面单元
  - 周期性边界条件
  - 已划分网格的装配体，可直接分析

### CSV导出（可选）
- **文件名**：`FiberCenters_VfXX_YYYYMMDD_HHMMSS.csv`
- **内容**：
  - 纤维中心坐标
  - RVE尺寸信息
  - 体积分数信息
  - 时间戳和元数据

**CSV格式示例**：
```csv
# RVE Fiber Center Coordinates
# Generated: 2025-10-20 14:30:00
# RVE Size (Width x Height): 1.000000 x 1.000000
# Fiber Radius: 0.050000
# Target Volume Fraction: 0.5000
# Total Fiber Count (effective in RVE): 127
Fiber_ID,X_Coordinate,Y_Coordinate
1,0.12345678,0.23456789
2,0.34567890,0.45678901
...
```

## 📖 算法详解

### 纤维布置策略

脚本采用复杂的三阶段算法：

1. **RSA播种阶段**：
   - 使用随机序列吸附生成初始纤维位置
   - 由`RSA_SEEDING_RATIO`参数控制
   - 确保初始布置无重叠

2. **锚定松弛阶段**：
   - 优化纤维位置以达到目标体积分数
   - 保持最小距离约束
   - 迭代调整位置实现最优填充

3. **强制校正阶段**：
   - 最终调整以满足精确体积分数要求
   - 处理边界条件和周期性约束

### 周期性边界条件

实现包括：
- 自动配对对边节点
- 位移连续性的方程约束
- 用于控制周期性变形的参考点
- 支持X和Y方向

### 内聚力界面实现

- 自动检测纤维-基体界面
- 插入COH2D4单元
- 支持混合模式断裂准则（BK准则）
- 粘性稳定化以提高数值稳定性

## 🎨 定制指南

### 调整体积分数
```python
TARGET_VF = 0.30  # 设置30%纤维体积分数
```
**注意**：可达到的最大Vf取决于纤维半径和RVE尺寸。如果无法满足目标，脚本会报错。

### 改变分布模式
```python
# 均匀分布，用于生产环境
RSA_SEEDING_RATIO = 0.9

# 真实团簇分布
RSA_SEEDING_RATIO = 0.1

# 平衡方案
RSA_SEEDING_RATIO = 0.5
```

### 修改网格密度
```python
# 更细网格（更精确，更慢）
GLOBAL_SEED_SIZE = 0.005

# 更粗网格（更快，精度较低）
GLOBAL_SEED_SIZE = 0.02
```

## ⚠️ 重要说明

### 限制条件
- 最大纤维体积分数受几何约束限制
- 内聚力单元插入可能增加计算成本
- 周期性边界条件需要谨慎施加载荷

### 最佳实践
1. 从较低体积分数(Vf < 0.6)开始测试
2. 根据网格密度需求调整`MIN_DIST_FACTOR`
3. 使用导出的CSV坐标验证纤维间距
4. 运行分析前检查网格质量
5. 为特定复合材料系统使用适当的材料参数

### 故障排除

**问题**：脚本报告"无法达到目标体积分数"
- **解决方案**：降低`TARGET_VF`或增加`RVE_SIZE`或减小`FIBER_RADIUS`

**问题**：网格生成失败
- **解决方案**：增加`GLOBAL_SEED_SIZE`或调整几何

**问题**：分析收敛问题
- **解决方案**：检查内聚力参数，调整`COHESIVE_STAB_COEFF`

## 📝 引用

如果您在研究中使用此代码，请引用：

```bibtex
@software{2d_rve_model_2025,
  author = {刘正鹏},
  title = {ABAQUS二维RVE建模工具},
  year = {2025},
  url = {https://github.com/ZPL-03/2D_RVE_Model},
  version = {1.0}
}
```

## 👥 作者

**刘正鹏 (Liu Zhengpeng)**
- GitHub: [@小盆i](https://github.com/ZPL-03)
- 邮箱：1370872708@qq.com / Zhengpeng0105@gmail.com
- 技术博客：CSDN/知乎 @小盆i

## 🤝 贡献

欢迎贡献！请阅读[CONTRIBUTING.md](CONTRIBUTING.md)了解行为准则和提交拉取请求的流程。

## 📄 许可证

本项目采用MIT许可证 - 详见[LICENSE](LICENSE)文件。

## 🙏 致谢

- 受复合材料微观力学研究社区的启发
- 感谢所有贡献者和用户的反馈和改进

## 📮 联系与支持

- **问题反馈**：通过[GitHub Issues](https://github.com/ZPL-03/2D_RVE_Model/issues)报告错误和请求功能
- **讨论交流**：在[GitHub Discussions](https://github.com/ZPL-03/2D_RVE_Model/discussions)参与讨论
- **私人咨询**：请发送邮件至 Zhengpeng0105@gmail.com

## 📸 示例展示

### 生成的RVE模型示例
（建议添加一些模型截图）

### 典型应用案例
（可以添加一些使用此工具发表的论文或项目）

## 🗺️ 发展路线图

- [ ] 支持3D RVE模型生成
- [ ] 添加椭圆形纤维支持
- [ ] 实现纤维取向控制
- [ ] 增加GUI界面
- [ ] 支持更多类型的边界条件
- [ ] 添加自动优化算法

## ⭐ Star History

如果这个项目对您有帮助，请给它一个星标⭐！

---

**版本**: v1.0  
**最后更新**: 2025-10-20  
**维护状态**: 活跃维护中
