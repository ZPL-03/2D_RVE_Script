# -*- coding: utf-8 -*-
# ##################################################################
#
#           可控制纤维体积分数的Abaqus二维RVE建模脚本
#
# 功能说明:
# 1. 根据用户设定的纤维体积分数(Vf)自动生成二维RVE几何模型
# 2. 采用混合算法（RSA播种 + 锚定松弛 + 强制校正）排布纤维，兼顾速度与分布质量
# 3. 保证纤维间的最小间距，无法满足时报错终止
# 4. 自动创建材料（基体+纤维+界面）、赋予截面、划分网格
# 5. 施加周期性边界条件(PBCs)实现XY方向周期性
# 6. 在纤维-基体界面插入二维cohesive单元(COH2D4)
# 7. 利用几何位置算法准确识别纤维和基体（包括边角碎片）
# 8. 导出纤维中心坐标为CSV文件，便于验证和后续分析
# 9. 验证纤维间距并输出详细统计信息
#
# 作者: 刘正鹏 (Liu Zhengpeng)
# 版本: v1.0
# 创建日期: 2025-09-30
# 最后更新: 2025-XX-XX
# 适用软件: ABAQUS 2023
# Python版本: 2.7 (ABAQUS内置)
# 技术交流: GitHub/CSDN/知乎 @小盆i
# 联系方式: 1370872708@qq.com / Zhengpeng0105@gmail.com
#
# ##################################################################

# ============ Abaqus相关模块导入 ============
from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup

# ============ Python标准库导入 ============
import math
import random as rd
import time
import os

# ============ Abaqus子模块导入 ============
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from sketch import *
from visualization import *
from connectorBehavior import *

executeOnCaeStartup()


# =================================================================
#                 CSV坐标输出模块
# =================================================================
def exportFiberCentersToCSV(fiber_centers, filename, rveSize, fiberRadius, target_Vf):
    """将纤维中心坐标导出为CSV文件

    参数:
        fiber_centers: 纤维中心坐标列表 [(x1,y1), (x2,y2), ...]
        filename: 输出的CSV文件名
        rveSize: RVE尺寸 [宽度, 高度]
        fiberRadius: 纤维半径
        target_Vf: 目标体积分数

    返回:
        filepath: 成功则返回文件完整路径，失败则返回None

    注意:
        导出的坐标为RVE内的有效纤维中心，不包括周期性重复的纤维
    """
    try:
        work_dir = os.getcwd()
        filepath = os.path.join(work_dir, filename)

        with open(filepath, 'w') as f:
            f.write("# RVE Fiber Center Coordinates\n")
            f.write("# Generated: %s\n" % time.strftime("%Y-%m-%d %H:%M:%S"))
            f.write("# RVE Size (Width x Height): %.6f x %.6f\n" % (rveSize[0], rveSize[1]))
            f.write("# Fiber Radius: %.6f\n" % fiberRadius)
            f.write("# Target Volume Fraction: %.4f\n" % target_Vf)
            f.write("# Total Fiber Count (effective in RVE): %d\n" % len(fiber_centers))
            f.write("# Note: These are original fiber centers, NOT including periodic images\n")
            f.write("#\n")
            f.write("Fiber_ID,X_Coordinate,Y_Coordinate\n")

            for i, (x, y) in enumerate(fiber_centers, start=1):
                f.write("%d,%.8f,%.8f\n" % (i, x, y))

        print("\n" + "=" * 60)
        print("SUCCESS: Fiber coordinates exported to CSV file")
        print("File location: %s" % filepath)
        print("Total fibers exported: %d" % len(fiber_centers))
        print("=" * 60 + "\n")
        return filepath

    except Exception as e:
        print("\nWARNING: Failed to export coordinates: %s\n" % str(e))
        return None


# =================================================================
#                 周期性边界条件 (PBCs) 辅助函数
# =================================================================
def createReferencePoints(model):
    """创建两个参考点用于施加周期性边界条件

    参数:
        model: Abaqus模型对象

    功能:
        创建Ref-LR和Ref-BT两个参考点，分别用于左右边界和上下边界的周期性约束
    """
    rootAssembly = model.rootAssembly

    p_LR = model.Part(name='Ref-LR-Part', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
    p_LR.ReferencePoint(point=(0.0, 0.0, 0.0))

    p_BT = model.Part(name='Ref-BT-Part', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
    p_BT.ReferencePoint(point=(0.0, 0.0, 0.0))

    rootAssembly.Instance(name='Ref-LR-Instance', part=p_LR, dependent=ON)
    rootAssembly.Instance(name='Ref-BT-Instance', part=p_BT, dependent=ON)

    rootAssembly.Set(name='set_RefPoint_LR',
                     referencePoints=(rootAssembly.instances['Ref-LR-Instance'].referencePoints[1],))
    rootAssembly.Set(name='set_RefPoint_BT',
                     referencePoints=(rootAssembly.instances['Ref-BT-Instance'].referencePoints[1],))


def getRVEDimensions(model, instanceName):
    """获取RVE模型的精确边界尺寸

    参数:
        model: Abaqus模型对象
        instanceName: 实例名称

    返回:
        (xMin, xMax, yMin, yMax): RVE的边界坐标
    """
    nodes = model.rootAssembly.instances[instanceName].nodes

    xMin = min([n.coordinates[0] for n in nodes])
    xMax = max([n.coordinates[0] for n in nodes])
    yMin = min([n.coordinates[1] for n in nodes])
    yMax = max([n.coordinates[1] for n in nodes])

    return xMin, xMax, yMin, yMax


def getBoundaryNodes(model, instanceName, dimensions):
    """获取并分类RVE模型四条边界上的所有节点

    参数:
        model: Abaqus模型对象
        instanceName: 实例名称
        dimensions: 边界尺寸 (xMin, xMax, yMin, yMax)

    返回:
        (nodes_left, nodes_right, nodes_bottom, nodes_top): 四条边界上的节点列表
    """
    nodes = model.rootAssembly.instances[instanceName].nodes
    nodes_left, nodes_right, nodes_bottom, nodes_top = [], [], [], []
    xMin, xMax, yMin, yMax = dimensions
    tol = 1e-6

    for n in nodes:
        x, y = n.coordinates[0], n.coordinates[1]
        if abs(x - xMin) < tol: nodes_left.append(n)
        if abs(x - xMax) < tol: nodes_right.append(n)
        if abs(y - yMin) < tol: nodes_bottom.append(n)
        if abs(y - yMax) < tol: nodes_top.append(n)

    nodes_left.sort(key=lambda n: n.coordinates[1])
    nodes_right.sort(key=lambda n: n.coordinates[1])
    nodes_bottom.sort(key=lambda n: n.coordinates[0])
    nodes_top.sort(key=lambda n: n.coordinates[0])

    return nodes_left, nodes_right, nodes_bottom, nodes_top


def pairBoundaryNodes(slave_nodes, master_nodes, tolerance, coord_index):
    """为两条相对边界上的节点进行配对

    参数:
        slave_nodes: 从属边界的节点列表
        master_nodes: 主边界的节点列表
        tolerance: 配对容差
        coord_index: 配对时使用的坐标索引（0为x，1为y）

    返回:
        paired_nodes: 配对结果列表 [(slave_node1, master_node1), ...]
    """
    paired_nodes = []
    master_pool = list(master_nodes)

    for s_node in slave_nodes:
        s_coord = s_node.coordinates
        best_match_node = None
        min_dist = float('inf')

        for m_node in master_pool:
            m_coord = m_node.coordinates
            if abs(s_coord[coord_index] - m_coord[coord_index]) < tolerance:
                dist = math.sqrt(sum([(sc - mc) ** 2 for sc, mc in zip(s_coord, m_coord)]))
                if dist < min_dist:
                    min_dist = dist
                    best_match_node = m_node

        if best_match_node:
            paired_nodes.append((s_node, best_match_node))
            master_pool.remove(best_match_node)

    return paired_nodes


def applyPeriodicConstraints(model, instanceName, node_pairs, pair_type):
    """为已配对的节点施加约束方程来实现周期性边界条件

    参数:
        model: Abaqus模型对象
        instanceName: 实例名称
        node_pairs: 节点配对列表
        pair_type: 配对类型（'Left-Right' 或 'Bottom-Top'）

    功能:
        为每对节点创建约束方程: u_slave - u_master + u_ref = 0
    """
    r_assy = model.rootAssembly
    inst = r_assy.instances[instanceName]

    if pair_type == 'Left-Right':
        ref_point_name = 'set_RefPoint_LR'
        tag1, tag2 = 'L', 'R'
        coeffs = (1.0, -1.0, 1.0)
    elif pair_type == 'Bottom-Top':
        ref_point_name = 'set_RefPoint_BT'
        tag1, tag2 = 'B', 'T'
        coeffs = (1.0, -1.0, 1.0)
    else:
        return

    for i, (node1, node2) in enumerate(node_pairs):
        set1_name = 'set_Node-%s-%d' % (tag1, i + 1)
        set2_name = 'set_Node-%s-%d' % (tag2, i + 1)
        r_assy.Set(nodes=inst.nodes.sequenceFromLabels(labels=(node1.label,)), name=set1_name)
        r_assy.Set(nodes=inst.nodes.sequenceFromLabels(labels=(node2.label,)), name=set2_name)

        model.Equation(name='Eq-%s%s-X-%d' % (tag1, tag2, i + 1),
                       terms=((coeffs[0], set1_name, 1),
                              (coeffs[1], set2_name, 1),
                              (coeffs[2], ref_point_name, 1)))

        model.Equation(name='Eq-%s%s-Y-%d' % (tag1, tag2, i + 1),
                       terms=((coeffs[0], set1_name, 2),
                              (coeffs[1], set2_name, 2),
                              (coeffs[2], ref_point_name, 2)))


# =================================================================
#                 纤维排布算法模块
# =================================================================
def _relax_coords_anchored(initial_coords, seeding_count, fiberCount, rveSize, fiberRadius, minDistance):
    """锚定松弛算法

    参数:
        initial_coords: 初始纤维坐标列表
        seeding_count: RSA播种的纤维数量（这些纤维作为锚点，移动受限）
        fiberCount: 总纤维数量
        rveSize: RVE尺寸 [宽度, 高度]
        fiberRadius: 纤维半径
        minDistance: 最小间距

    返回:
        coords: 松弛后的纤维坐标列表

    功能:
        通过模拟斥力来调整纤维位置，锚点纤维的移动受到阻尼限制
    """
    print("--- Initializing Anchored Relaxation Process ---")

    coords = [list(c) for c in initial_coords]
    max_iterations = 2000
    movement_factor = 0.5
    anchor_damping_factor = 0.05
    min_dist_sq = minDistance ** 2

    for iter_num in range(max_iterations):
        max_movement_sq = 0.0
        net_forces = [[0.0, 0.0] for _ in range(fiberCount)]

        for i in range(fiberCount):
            for j in range(i + 1, fiberCount):
                dx = coords[j][0] - coords[i][0]
                dy = coords[j][1] - coords[i][1]

                if dx > rveSize[0] / 2: dx -= rveSize[0]
                if dx < -rveSize[0] / 2: dx += rveSize[0]
                if dy > rveSize[1] / 2: dy -= rveSize[1]
                if dy < -rveSize[1] / 2: dy += rveSize[1]

                dist_sq = dx * dx + dy * dy

                if dist_sq < min_dist_sq:
                    dist = math.sqrt(dist_sq) if dist_sq > 0 else 1e-9
                    overlap = minDistance - dist
                    force_magnitude = overlap

                    force_x = force_magnitude * (dx / dist)
                    force_y = force_magnitude * (dy / dist)

                    net_forces[i][0] -= force_x
                    net_forces[i][1] -= force_y
                    net_forces[j][0] += force_x
                    net_forces[j][1] += force_y

        if not any(any(f) for f in net_forces):
            print("--- No overlaps detected. System stable. ---")
            return coords

        for i in range(fiberCount):
            move_x = net_forces[i][0] * movement_factor
            move_y = net_forces[i][1] * movement_factor

            if i < seeding_count:
                move_x *= anchor_damping_factor
                move_y *= anchor_damping_factor

            coords[i][0] += move_x
            coords[i][1] += move_y

            coords[i][0] %= rveSize[0]
            coords[i][1] %= rveSize[1]

            current_movement_sq = move_x ** 2 + move_y ** 2
            if current_movement_sq > max_movement_sq:
                max_movement_sq = current_movement_sq

        if (iter_num + 1) % 50 == 0:
            print("... Relaxation Iteration %d, Max movement: %.2e" %
                  (iter_num + 1, math.sqrt(max_movement_sq)))

        if max_movement_sq < (1e-6 * fiberRadius) ** 2:
            print("--- Converged after %d iterations. ---" % (iter_num + 1))
            return coords

    print("--- Max relaxation iterations reached. ---")
    return coords


def _final_check_and_enforce(coords_in, fiberCount, rveSize, minDistance):
    """最终检查和强制校正

    参数:
        coords_in: 输入的纤维坐标列表
        fiberCount: 纤维数量
        rveSize: RVE尺寸
        minDistance: 最小间距

    返回:
        coords: 校正后的坐标列表

    功能:
        检查所有纤维对的间距，对违反最小间距要求的纤维进行强制移动
        如果无法满足最小间距要求，抛出异常
    """
    coords = [list(c) for c in coords_in]
    min_dist_sq = minDistance ** 2
    max_correction_iter = 50000

    for iter_num in range(max_correction_iter):
        min_dist_sq_found = float('inf')
        worst_offenders = None

        for i in range(fiberCount):
            for j in range(i + 1, fiberCount):
                dx = coords[j][0] - coords[i][0]
                dy = coords[j][1] - coords[i][1]

                if dx > rveSize[0] / 2: dx -= rveSize[0]
                if dx < -rveSize[0] / 2: dx += rveSize[0]
                if dy > rveSize[1] / 2: dy -= rveSize[1]
                if dy < -rveSize[1] / 2: dy += rveSize[1]

                dist_sq = dx * dx + dy * dy

                if dist_sq < min_dist_sq_found:
                    min_dist_sq_found = dist_sq
                    worst_offenders = (i, j, dx, dy)

        if min_dist_sq_found >= min_dist_sq:
            print("--- Final Check PASSED after %d iterations. ---" % (iter_num + 1))
            return coords

        if worst_offenders:
            i, j, dx, dy = worst_offenders
            dist = math.sqrt(min_dist_sq_found) if min_dist_sq_found > 0 else 1e-9
            overlap = minDistance - dist
            move_dist = overlap / 2.0 + 1e-8

            move_x = move_dist * (dx / dist)
            move_y = move_dist * (dy / dist)

            coords[i][0] -= move_x
            coords[i][1] -= move_y
            coords[j][0] += move_x
            coords[j][1] += move_y

            coords[i][0] %= rveSize[0]
            coords[i][1] %= rveSize[1]
            coords[j][0] %= rveSize[0]
            coords[j][1] %= rveSize[1]

        if (iter_num + 1) % 100 == 0:
            print("... Correction Iteration %d, Min distance: %.6f (Target: %.6f)" %
                  (iter_num + 1, math.sqrt(min_dist_sq_found), minDistance))

    final_min_dist_sq = float('inf')
    violating_pairs = 0
    violating_fibers = set()

    for i in range(fiberCount):
        for j in range(i + 1, fiberCount):
            dx = coords[j][0] - coords[i][0]
            dy = coords[j][1] - coords[i][1]
            if dx > rveSize[0] / 2: dx -= rveSize[0]
            if dx < -rveSize[0] / 2: dx += rveSize[0]
            if dy > rveSize[1] / 2: dy -= rveSize[1]
            if dy < -rveSize[1] / 2: dy += rveSize[1]
            dist_sq = dx * dx + dy * dy

            if dist_sq < final_min_dist_sq:
                final_min_dist_sq = dist_sq

            if dist_sq < min_dist_sq:
                violating_pairs += 1
                violating_fibers.add(i)
                violating_fibers.add(j)

    if final_min_dist_sq < min_dist_sq:
        error_msg = (
                "FATAL ERROR: Cannot satisfy minimum distance after %d iterations.\n"
                " -> Final min distance: %.6f (Target: >= %.6f)\n"
                " -> Violating pairs: %d, Involved fibers: %d\n"
                " -> Please reduce Vf or increase RVE size."
                % (max_correction_iter, math.sqrt(final_min_dist_sq), minDistance,
                   violating_pairs, len(violating_fibers))
        )
        raise Exception(error_msg)
    else:
        print("--- Final Check PASSED. ---")
        return coords


# =================================================================
#         纤维-基体分类算法
# =================================================================
def buildAllFiberCenters(fiber_centers, rveSize, fiberRadius):
    """构建包含周期性镜像的完整纤维中心列表

    参数:
        fiber_centers: 原始RVE内的纤维中心列表 [(x1,y1), (x2,y2), ...]
        rveSize: RVE尺寸 [width, height]
        fiberRadius: 纤维半径

    返回:
        unique_centers: 包含所有周期性镜像的纤维中心列表

    功能:
        为靠近边界的纤维创建周期性镜像，用于准确判断面的归属
    """
    all_fiber_centers = []

    for xt, yt in fiber_centers:
        all_fiber_centers.append((xt, yt))

        if xt < fiberRadius:
            all_fiber_centers.append((xt + rveSize[0], yt))
        if xt > rveSize[0] - fiberRadius:
            all_fiber_centers.append((xt - rveSize[0], yt))

        if yt < fiberRadius:
            all_fiber_centers.append((xt, yt + rveSize[1]))
        if yt > rveSize[1] - fiberRadius:
            all_fiber_centers.append((xt, yt - rveSize[1]))

        if xt < fiberRadius and yt < fiberRadius:
            dist_to_corner = math.sqrt(xt ** 2 + yt ** 2)
            if dist_to_corner < fiberRadius:
                all_fiber_centers.append((xt + rveSize[0], yt + rveSize[1]))

        if xt < fiberRadius and yt > rveSize[1] - fiberRadius:
            dist_to_corner = math.sqrt(xt ** 2 + (rveSize[1] - yt) ** 2)
            if dist_to_corner < fiberRadius:
                all_fiber_centers.append((xt + rveSize[0], yt - rveSize[1]))

        if xt > rveSize[0] - fiberRadius and yt > rveSize[1] - fiberRadius:
            dist_to_corner = math.sqrt((rveSize[0] - xt) ** 2 + (rveSize[1] - yt) ** 2)
            if dist_to_corner < fiberRadius:
                all_fiber_centers.append((xt - rveSize[0], yt - rveSize[1]))

        if xt > rveSize[0] - fiberRadius and yt < fiberRadius:
            dist_to_corner = math.sqrt((rveSize[0] - xt) ** 2 + yt ** 2)
            if dist_to_corner < fiberRadius:
                all_fiber_centers.append((xt - rveSize[0], yt + rveSize[1]))

    unique_centers = []
    tolerance = 1e-6
    for center in all_fiber_centers:
        is_duplicate = False
        for existing in unique_centers:
            if abs(center[0] - existing[0]) < tolerance and abs(center[1] - existing[1]) < tolerance:
                is_duplicate = True
                break
        if not is_duplicate:
            unique_centers.append(center)

    return unique_centers


def getFaceCenterFromVertices(face):
    """通过面的顶点坐标计算几何中心

    参数:
        face: Abaqus面对象

    返回:
        (x, y): 面的几何中心坐标，如果失败则返回None
    """
    try:
        vertices = face.getVertices()
        if not vertices or len(vertices) == 0:
            return None

        x_coords = []
        y_coords = []

        for vertex in vertices:
            try:
                if hasattr(vertex, 'pointOn'):
                    coord = vertex.pointOn[0]
                    x_coords.append(coord[0])
                    y_coords.append(coord[1])
            except:
                continue

        if len(x_coords) > 0 and len(y_coords) > 0:
            center_x = sum(x_coords) / len(x_coords)
            center_y = sum(y_coords) / len(y_coords)
            return (center_x, center_y)
        else:
            return None

    except Exception as e:
        print("      Debug: Error in getFaceCenterFromVertices: %s" % str(e))
        return None


def getFaceCenterAlternative(face):
    """通过面的pointOn属性获取中心点

    参数:
        face: Abaqus面对象

    返回:
        (x, y): 面上的一个点坐标，如果失败则返回None
    """
    try:
        if hasattr(face, 'pointOn'):
            point_on_face = face.pointOn
            if point_on_face and len(point_on_face) > 0:
                coord = point_on_face[0]
                if len(coord) >= 2:
                    return (coord[0], coord[1])
        return None
    except:
        return None


def classifyFacesImproved(all_faces, fiber_centers, rveSize, fiberRadius, rveArea):
    """基于几何位置的面分类算法

    参数:
        all_faces: 所有面的列表
        fiber_centers: 原始纤维中心坐标列表
        rveSize: RVE尺寸
        fiberRadius: 纤维半径
        rveArea: RVE面积

    返回:
        (fiber_faces_list, matrix_faces_list): 纤维面列表和基体面列表

    功能:
        通过计算面中心到纤维中心的距离来判断面的类型
        如果距离小于纤维半径，则为纤维面；否则为基体面
    """
    print("\n  === Starting Improved Face Classification ===")
    print("  Total faces to classify: %d" % len(all_faces))
    print("  Original fiber count: %d" % len(fiber_centers))

    sorted_faces = sorted(all_faces, key=lambda f: f.getSize(), reverse=True)
    print("\n  Step 1: Area-based sorting")
    print("    Largest face area: %.6e (assumed to be main matrix)" % sorted_faces[0].getSize())

    all_fiber_centers = buildAllFiberCenters(fiber_centers, rveSize, fiberRadius)
    print("\n  Step 2: Building complete fiber center list")
    print("    Original fibers: %d" % len(fiber_centers))
    print("    Total fiber centers (with periodicity): %d" % len(all_fiber_centers))

    matrix_faces_list = [sorted_faces[0]]
    potential_faces = sorted_faces[1:]

    print("\n  Step 3: Initial classification")
    print("    Main matrix face: 1")
    print("    Faces to classify: %d" % len(potential_faces))

    print("\n  Step 4: Detailed classification using geometry")

    if len(potential_faces) < len(all_fiber_centers):
        error_msg = (
                "FATAL ERROR: Face count anomaly detected!\n"
                " -> Potential faces: %d\n"
                " -> all fiber count: %d\n"
                % (len(potential_faces), len(all_fiber_centers))
        )
        raise Exception(error_msg)

    elif len(potential_faces) == len(all_fiber_centers):
        print("    Potential face count equals all fiber count (%d = %d)" %
              (len(potential_faces), len(all_fiber_centers)))
        print("    Perfect match - all potential faces are fibers")
        fiber_faces_list = potential_faces
        matrix_fragment_count = 0

    else:
        print("    Potential face count > all fiber count (%d > %d)" %
              (len(potential_faces), len(all_fiber_centers)))
        print("    Matrix fragments detected - performing geometry-based classification")

        fiber_faces_list = []
        matrix_fragment_count = 0

        for idx, face in enumerate(potential_faces):
            face_center = None

            face_center = getFaceCenterAlternative(face)

            if face_center is None:
                try:
                    centroid = face.getCentroid()
                    if centroid is not None:
                        if len(centroid) >= 2:
                            face_center = (centroid[0], centroid[1])
                        elif len(centroid) == 3:
                            face_center = (centroid[0], centroid[1])
                except:
                    pass

            if face_center is None:
                face_center = getFaceCenterFromVertices(face)

            if face_center is None:
                error_msg = (
                        "FATAL ERROR: Unable to determine face center for face %d\n"
                        " -> All center calculation methods failed\n"
                        " -> Face area: %.6e\n"
                        " -> Cannot proceed with classification"
                        % (idx, face.getSize())
                )
                raise Exception(error_msg)

            face_x, face_y = face_center[0], face_center[1]

            min_dist = float('inf')
            for fc_x, fc_y in all_fiber_centers:
                dist = math.sqrt((face_x - fc_x) ** 2 + (face_y - fc_y) ** 2)
                if dist < min_dist:
                    min_dist = dist

            if min_dist < fiberRadius:
                fiber_faces_list.append(face)
            else:
                matrix_faces_list.append(face)
                matrix_fragment_count += 1
                print("      [Fragment %d] Area=%.2e, CenterDist=%.6f > Radius=%.6f" %
                      (matrix_fragment_count, face.getSize(), min_dist, fiberRadius))

    print("\n  === Classification Complete ===")
    print("  Fiber faces identified: %d" % len(fiber_faces_list))
    print("  Matrix faces total: %d" % len(matrix_faces_list))

    fiber_total_area = sum([f.getSize() for f in fiber_faces_list])
    actual_Vf = fiber_total_area / rveArea
    target_Vf_from_count = (len(fiber_centers) * math.pi * fiberRadius ** 2) / rveArea

    print("\n  Validation:")
    print("    Target Vf (from fiber count): %.4f" % target_Vf_from_count)
    print("    Actual Vf (from classified faces): %.4f" % actual_Vf)
    print("    Deviation: %.2f%%" % (abs(actual_Vf - target_Vf_from_count) / target_Vf_from_count * 100))

    return fiber_faces_list, matrix_faces_list


# =================================================================
#                 纤维间距验证模块
# =================================================================
def verifyMinimumFiberDistance(fiber_centers, rveSize, fiberRadius, minDistanceFactor):
    """验证所有纤维对之间的距离是否满足最小间距要求

    参数:
        fiber_centers: 纤维中心坐标列表 [(x1,y1), (x2,y2), ...]
        rveSize: RVE尺寸 [宽度, 高度]
        fiberRadius: 纤维半径
        minDistanceFactor: 最小间距因子（例如2.05）

    返回:
        verification_passed: 布尔值，是否通过验证
        stats: 统计信息字典
    """
    print("\n" + "=" * 70)
    print("FIBER DISTANCE VERIFICATION")
    print("=" * 70)

    if len(fiber_centers) == 0:
        print("No fibers to verify.")
        return True, {}

    if len(fiber_centers) == 1:
        print("Only 1 fiber - no distance check needed.")
        return True, {'fiber_count': 1}

    fiberCount = len(fiber_centers)
    minDistance = minDistanceFactor * fiberRadius

    print("\nConfiguration:")
    print("  Total Fibers: %d" % fiberCount)
    print("  Fiber Radius: %.6f mm" % fiberRadius)
    print("  Min Distance Factor: %.2f" % minDistanceFactor)
    print("  Required Min Distance: %.6f mm (%.2f x %.6f)" %
          (minDistance, minDistanceFactor, fiberRadius))
    print("  RVE Size: %.6f x %.6f mm" % (rveSize[0], rveSize[1]))

    print("\nCalculating inter-fiber distances...")

    min_distance_found = float('inf')
    max_distance_found = 0.0
    total_pairs = 0
    violating_pairs = []
    all_distances = []

    for i in range(fiberCount):
        for j in range(i + 1, fiberCount):
            dx = fiber_centers[j][0] - fiber_centers[i][0]
            dy = fiber_centers[j][1] - fiber_centers[i][1]

            if dx > rveSize[0] / 2.0:
                dx -= rveSize[0]
            if dx < -rveSize[0] / 2.0:
                dx += rveSize[0]
            if dy > rveSize[1] / 2.0:
                dy -= rveSize[1]
            if dy < -rveSize[1] / 2.0:
                dy += rveSize[1]

            distance = math.sqrt(dx * dx + dy * dy)
            all_distances.append(distance)
            total_pairs += 1

            if distance < min_distance_found:
                min_distance_found = distance
            if distance > max_distance_found:
                max_distance_found = distance

            if distance < minDistance:
                violating_pairs.append({
                    'fiber_i': i + 1,
                    'fiber_j': j + 1,
                    'distance': distance,
                    'violation': minDistance - distance,
                    'center_i': fiber_centers[i],
                    'center_j': fiber_centers[j]
                })

    avg_distance = sum(all_distances) / len(all_distances)
    all_distances.sort()
    median_distance = all_distances[len(all_distances) // 2]

    print("\n" + "-" * 70)
    print("DISTANCE STATISTICS")
    print("-" * 70)
    print("  Total Fiber Pairs Checked: %d" % total_pairs)
    print("  Minimum Distance Found: %.6f mm" % min_distance_found)
    print("  Maximum Distance Found: %.6f mm" % max_distance_found)
    print("  Average Distance: %.6f mm" % avg_distance)
    print("  Median Distance: %.6f mm" % median_distance)

    distance_ratio = min_distance_found / minDistance
    print("\n  Distance Ratio (min_found / required): %.4f" % distance_ratio)

    if distance_ratio >= 1.0:
        print("  Status: PASSED (%.2f%% above requirement)" %
              ((distance_ratio - 1.0) * 100))
    else:
        print("  Status: FAILED (%.2f%% below requirement)" %
              ((1.0 - distance_ratio) * 100))

    print("\n" + "-" * 70)
    if len(violating_pairs) == 0:
        print("VERIFICATION RESULT: ALL DISTANCES SATISFY MINIMUM REQUIREMENT")
        print("-" * 70)
        verification_passed = True
    else:
        print("VERIFICATION RESULT: MINIMUM DISTANCE VIOLATIONS DETECTED")
        print("-" * 70)
        print("  Number of Violating Pairs: %d" % len(violating_pairs))
        print("\nTop 10 Violations (sorted by severity):")
        print("-" * 70)

        violating_pairs.sort(key=lambda x: x['violation'], reverse=True)

        for idx, violation in enumerate(violating_pairs[:10], 1):
            print("\n  Violation #%d:" % idx)
            print("    Fiber Pair: #%d ↔ #%d" %
                  (violation['fiber_i'], violation['fiber_j']))
            print("    Actual Distance: %.6f mm" % violation['distance'])
            print("    Required Distance: %.6f mm" % minDistance)
            print("    Shortfall: %.6f mm (%.2f%% below requirement)" %
                  (violation['violation'],
                   violation['violation'] / minDistance * 100))
            print("    Center #%d: (%.6f, %.6f)" %
                  (violation['fiber_i'],
                   violation['center_i'][0],
                   violation['center_i'][1]))
            print("    Center #%d: (%.6f, %.6f)" %
                  (violation['fiber_j'],
                   violation['center_j'][0],
                   violation['center_j'][1]))

        if len(violating_pairs) > 10:
            print("\n  ... and %d more violations" %
                  (len(violating_pairs) - 10))

        verification_passed = False

    print("=" * 70 + "\n")

    stats = {
        'fiber_count': fiberCount,
        'total_pairs': total_pairs,
        'min_distance': min_distance_found,
        'max_distance': max_distance_found,
        'avg_distance': avg_distance,
        'median_distance': median_distance,
        'required_distance': minDistance,
        'distance_ratio': distance_ratio,
        'violations_count': len(violating_pairs),
        'violations': violating_pairs
    }

    return verification_passed, stats


# =================================================================
#                 主建模函数
# =================================================================
def createRVEModel(modelName='Model-1',
                   rveSize=[0.057, 0.057],
                   fiberRadius=0.0035,
                   target_Vf=0.5,
                   minDistanceFactor=2.05,
                   globalSeedSize=0.0005,
                   pairingToleranceFactor=0.1,
                   rsa_seeding_ratio=0.9,
                   export_coordinates=True,
                   csv_filename=None,
                   fiber_E=15000.0,
                   fiber_nu=0.3,
                   matrix_E=3170.0,
                   matrix_nu=0.35,
                   matrix_friction_angle=16.0,
                   matrix_flow_stress_ratio=1.0,
                   matrix_dilation_angle=16.0,
                   matrix_hardening_yield=106.4,
                   matrix_hardening_plastic_strain=0.0,
                   matrix_damage_strain=0.01,
                   matrix_damage_stress_triax=0.0,
                   matrix_damage_strain_rate=0.0,
                   matrix_damage_displacement=5e-05,
                   cohesive_K_nn=1e8,
                   cohesive_K_ss=1e8,
                   cohesive_K_tt=1e8,
                   cohesive_t_n=44.0,
                   cohesive_t_s=82.0,
                   cohesive_t_t=82.0,
                   cohesive_GIC=0.001,
                   cohesive_GIIC=0.002,
                   cohesive_GIIIC=0.002,
                   cohesive_eta=1.5,
                   cohesive_stab_coeff=0.0001):
    """创建RVE模型的主函数

    参数说明:
        modelName: 模型名称
        rveSize: RVE尺寸 [宽度, 高度]
        fiberRadius: 纤维半径
        target_Vf: 目标体积分数
        minDistanceFactor: 最小间距因子（实际最小间距 = 因子 × 纤维半径）
        globalSeedSize: 全局网格尺寸
        pairingToleranceFactor: 边界节点配对容差因子
        rsa_seeding_ratio: RSA播种比例（0-1，值越大分布越均匀）
        export_coordinates: 是否导出坐标到CSV文件
        csv_filename: CSV文件名（None则自动生成）

        纤维材料参数 (MPa):
        fiber_E: 弹性模量
        fiber_nu: 泊松比

        基体材料参数 (MPa):
        matrix_E: 弹性模量
        matrix_nu: 泊松比
        matrix_friction_angle: Drucker-Prager摩擦角
        matrix_flow_stress_ratio: Drucker-Prager流动应力比
        matrix_dilation_angle: Drucker-Prager膨胀角
        matrix_hardening_yield: 硬化屈服应力
        matrix_hardening_plastic_strain: 硬化对应塑性应变
        matrix_damage_strain: 韧性损伤起始应变
        matrix_damage_stress_triax: 应力三轴度
        matrix_damage_strain_rate: 应变率
        matrix_damage_displacement: 损伤演化位移

        界面材料参数:
        cohesive_K_nn, cohesive_K_ss, cohesive_K_tt: 界面刚度 (N/mm^3)
        cohesive_t_n, cohesive_t_s, cohesive_t_t: 界面强度 (MPa)
        cohesive_GIC, cohesive_GIIC, cohesive_GIIIC: 断裂能 (N/mm)
        cohesive_eta: BK准则指数
        cohesive_stab_coeff: 粘性稳定系数
    """

    # ==================== 步骤 1: 计算参数和纤维坐标 ====================
    print("Step 1: Calculating parameters and generating fiber coordinates...")

    rveArea = rveSize[0] * rveSize[1]
    fiberArea = math.pi * fiberRadius ** 2
    fiberCount = int(round((target_Vf * rveArea) / fiberArea))
    minDistance = minDistanceFactor * fiberRadius
    rsa_max_attempts = 500

    print("Target Vf: %.4f" % target_Vf)
    print("Calculated Fiber Count: %d" % fiberCount)
    print("RSA Seeding Ratio: %.2f" % rsa_seeding_ratio)
    print("Min Distance: %.6f (%.2fx radius)" % (minDistance, minDistanceFactor))

    if fiberCount == 0:
        print("Warning: Fiber count is zero.")
        fiber_centers = []
    elif fiberArea >= rveArea:
        print("Error: Single fiber area >= RVE area. Aborting.")
        return
    else:
        print("\n--- Stage 1: RSA Seeding ---")
        seeding_count = int(fiberCount * rsa_seeding_ratio)
        seeded_coords = []

        for i in range(seeding_count):
            placed = False
            for _ in range(rsa_max_attempts):
                xt = rd.uniform(0, rveSize[0])
                yt = rd.uniform(0, rveSize[1])
                is_too_close = False

                for xc, yc in seeded_coords:
                    dx = abs(xt - xc)
                    dy = abs(yt - yc)
                    if dx > rveSize[0] / 2: dx = rveSize[0] - dx
                    if dy > rveSize[1] / 2: dy = rveSize[1] - dy

                    if dx * dx + dy * dy < minDistance ** 2:
                        is_too_close = True
                        break

                if not is_too_close:
                    seeded_coords.append((xt, yt))
                    placed = True
                    break

            if not placed:
                print("RSA congested at %d fibers" % len(seeded_coords))
                break

        seeding_count_actual = len(seeded_coords)
        print("RSA placed %d anchor fibers" % seeding_count_actual)

        print("\n--- Stage 2: Random Placement ---")
        remaining = fiberCount - seeding_count_actual
        print("Placing %d fluid fibers..." % remaining)
        initial_coords = list(seeded_coords)

        for _ in range(remaining):
            initial_coords.append((rd.uniform(0, rveSize[0]),
                                   rd.uniform(0, rveSize[1])))

        print("\n--- Stage 3: Anchored Relaxation ---")
        relaxed_coords = _relax_coords_anchored(
            initial_coords, seeding_count_actual, fiberCount,
            rveSize, fiberRadius, minDistance
        )

        print("\n--- Stage 4: Final Verification ---")
        try:
            fiber_centers = _final_check_and_enforce(
                relaxed_coords, fiberCount, rveSize, minDistance
            )
        except Exception as e:
            print("\n" + "#" * 70)
            print(str(e))
            print("#" * 70 + "\n")
            return

    if export_coordinates and fiber_centers:
        if csv_filename is None:
            csv_filename = "FiberCenters_Vf%d_%s.csv" % (
                int(target_Vf * 100), time.strftime("%Y%m%d_%H%M%S")
            )
        exportFiberCentersToCSV(fiber_centers, csv_filename,
                                rveSize, fiberRadius, target_Vf)

    if fiber_centers:
        verification_passed, verification_stats = verifyMinimumFiberDistance(
            fiber_centers, rveSize, fiberRadius, minDistanceFactor
        )

        if not verification_passed:
            print("\nWARNING: Distance verification failed!")
            raise Exception("Minimum distance requirement not satisfied!")

    # ==================== 步骤 2: 创建几何 ====================
    print("\nStep 2: Creating RVE geometry...")

    xCoords, yCoords = [], []

    for xt, yt in fiber_centers:
        points_to_add = [(xt, yt)]

        if xt < fiberRadius:
            points_to_add.append((xt + rveSize[0], yt))
        if xt > rveSize[0] - fiberRadius:
            points_to_add.append((xt - rveSize[0], yt))
        if yt < fiberRadius:
            points_to_add.append((xt, yt + rveSize[1]))
        if yt > rveSize[1] - fiberRadius:
            points_to_add.append((xt, yt - rveSize[1]))
        if xt < fiberRadius and yt > rveSize[1] - fiberRadius:
            points_to_add.append((xt + rveSize[0], yt - rveSize[1]))
        if xt < fiberRadius and yt < fiberRadius:
            points_to_add.append((xt + rveSize[0], yt + rveSize[1]))
        if xt > rveSize[0] - fiberRadius and yt > rveSize[1] - fiberRadius:
            points_to_add.append((xt - rveSize[0], yt - rveSize[1]))
        if xt > rveSize[0] - fiberRadius and yt < fiberRadius:
            points_to_add.append((xt - rveSize[0], yt + rveSize[1]))

        unique_points = sorted(list(set(points_to_add)))
        for p in unique_points:
            xCoords.append(p[0])
            yCoords.append(p[1])

    if modelName in mdb.models:
        del mdb.models[modelName]

    mdb.Model(name=modelName, modelType=STANDARD_EXPLICIT)
    a = mdb.models[modelName]

    if 'Model-1' in mdb.models and modelName != 'Model-1':
        if len(mdb.models['Model-1'].parts) == 0:
            del mdb.models['Model-1']

    s_matrix = a.ConstrainedSketch(name='MatrixSketch', sheetSize=max(rveSize))
    s_matrix.rectangle(point1=(0.0, 0.0), point2=(rveSize[0], rveSize[1]))
    rvePart = a.Part(name='RVE-2D', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
    rvePart.BaseShell(sketch=s_matrix)

    if fiberCount > 0:
        s_partition = a.ConstrainedSketch(name='PartitionSketch', sheetSize=max(rveSize) * 3)

        for i in range(len(xCoords)):
            s_partition.CircleByCenterPerimeter(
                center=(xCoords[i], yCoords[i]),
                point1=(xCoords[i] + fiberRadius, yCoords[i])
            )

        face_to_partition = rvePart.faces[0]
        rvePart.PartitionFaceBySketch(faces=face_to_partition, sketch=s_partition)

    print("Step 2 Complete.")

    # ==================== 步骤 3: 创建集合 ====================
    print("\nStep 3: Creating geometry sets...")

    all_faces = rvePart.faces
    rvePart.Set(faces=all_faces, name='set_AllFaces')

    if len(all_faces) > 1 and fiber_centers:
        fiber_faces_list, matrix_faces_list = classifyFacesImproved(
            all_faces, fiber_centers, rveSize, fiberRadius, rveArea
        )

        if matrix_faces_list:
            rvePart.Set(name='set_MatrixFaces', faces=FaceArray(matrix_faces_list))
        else:
            rvePart.Set(name='set_MatrixFaces', faces=FaceArray())

        if fiber_faces_list:
            rvePart.Set(name='set_FiberFaces', faces=FaceArray(fiber_faces_list))
        else:
            rvePart.Set(name='set_FiberFaces', faces=FaceArray())
    else:
        rvePart.Set(name='set_MatrixFaces',
                    faces=all_faces if len(all_faces) == 1 else FaceArray())
        rvePart.Set(name='set_FiberFaces', faces=FaceArray())

    all_edges = rvePart.edges
    rvePart.Set(edges=all_edges, name='set_AllEdges')

    edge_left = all_edges.getByBoundingBox(0, 0, 0, 0, rveSize[1], 0)
    edge_right = all_edges.getByBoundingBox(rveSize[0], 0, 0, rveSize[0], rveSize[1], 0)
    edge_top = all_edges.getByBoundingBox(0, rveSize[1], 0, rveSize[0], rveSize[1], 0)
    edge_bottom = all_edges.getByBoundingBox(0, 0, 0, rveSize[0], 0, 0)

    rvePart.Set(edges=edge_left + edge_right + edge_top + edge_bottom, name='set_OuterEdges')

    rvePart.SetByBoolean(name='set_CohesiveEdges',
                         sets=(rvePart.sets['set_AllEdges'], rvePart.sets['set_OuterEdges']),
                         operation=DIFFERENCE)

    print("Step 3 Complete.")

    # ==================== 步骤 4: 材料和截面 ====================
    print("\nStep 4: Defining materials and sections...")

    matrixMaterial = a.Material(name='Material-Matrix')
    matrixMaterial.Elastic(table=((matrix_E, matrix_nu),))
    matrixMaterial.DruckerPrager(table=((matrix_friction_angle,
                                         matrix_flow_stress_ratio,
                                         matrix_dilation_angle),))
    matrixMaterial.druckerPrager.DruckerPragerHardening(
        table=((matrix_hardening_yield, matrix_hardening_plastic_strain),))
    matrixMaterial.DuctileDamageInitiation(
        table=((matrix_damage_strain,
                matrix_damage_stress_triax,
                matrix_damage_strain_rate),))
    matrixMaterial.ductileDamageInitiation.DamageEvolution(
        type=DISPLACEMENT, table=((matrix_damage_displacement,),))

    fiberMaterial = a.Material(name='Material-Fiber')
    fiberMaterial.Elastic(table=((fiber_E, fiber_nu),))

    cohesiveMaterial = a.Material(name='Material-Cohesive')
    cohesiveMaterial.Elastic(type=TRACTION,
                             table=((cohesive_K_nn, cohesive_K_ss, cohesive_K_tt),))
    cohesiveMaterial.QuadsDamageInitiation(
        table=((cohesive_t_n, cohesive_t_s, cohesive_t_t),))
    cohesiveMaterial.quadsDamageInitiation.DamageEvolution(
        type=ENERGY, mixedModeBehavior=BK, power=cohesive_eta,
        table=((cohesive_GIC, cohesive_GIIC, cohesive_GIIIC),))
    cohesiveMaterial.quadsDamageInitiation.DamageStabilizationCohesive(
        cohesiveCoeff=cohesive_stab_coeff)

    a.HomogeneousSolidSection(name='Section-Fiber',
                              material='Material-Fiber', thickness=None)
    a.HomogeneousSolidSection(name='Section-Matrix',
                              material='Material-Matrix', thickness=None)
    a.CohesiveSection(name='Section-Cohesive',
                      material='Material-Cohesive',
                      response=TRACTION_SEPARATION,
                      outOfPlaneThickness=None)

    if 'set_FiberFaces' in rvePart.sets and rvePart.sets['set_FiberFaces'].faces:
        rvePart.SectionAssignment(region=rvePart.sets['set_FiberFaces'],
                                  sectionName='Section-Fiber', offset=0.0)

    if 'set_MatrixFaces' in rvePart.sets and rvePart.sets['set_MatrixFaces'].faces:
        rvePart.SectionAssignment(region=rvePart.sets['set_MatrixFaces'],
                                  sectionName='Section-Matrix', offset=0.0)

    print("Step 4 Complete.")

    # ==================== 步骤 5: 网格 ====================
    print("\nStep 5: Meshing RVE...")

    rvePart.seedPart(size=globalSeedSize, deviationFactor=0.1, minSizeFactor=0.1)

    elemType_bulk = ElemType(elemCode=CPE4, elemLibrary=STANDARD,
                             elemDeletion=ON, viscosity=0.0001, maxDegradation=0.99)
    elemType_bulk_tri = ElemType(elemCode=CPE3, elemLibrary=STANDARD)

    rvePart.setElementType(regions=(rvePart.sets['set_AllFaces'].faces,),
                           elemTypes=(elemType_bulk, elemType_bulk_tri))

    rvePart.generateMesh()

    rvePart.insertElements(edges=rvePart.sets['set_CohesiveEdges'])

    all_elements = rvePart.elements
    rvePart.Set(elements=all_elements, name='set_AllElements')

    if 'set_FiberFaces' in rvePart.sets and rvePart.sets['set_FiberFaces'].faces:
        rvePart.Set(elements=rvePart.sets['set_FiberFaces'].elements,
                    name='set_FiberElements')
    else:
        rvePart.Set(elements=ElementArray([]), name='set_FiberElements')

    if 'set_MatrixFaces' in rvePart.sets and rvePart.sets['set_MatrixFaces'].faces:
        rvePart.Set(elements=rvePart.sets['set_MatrixFaces'].elements,
                    name='set_MatrixElements')
    else:
        rvePart.Set(elements=ElementArray([]), name='set_MatrixElements')

    rvePart.SetByBoolean(name='set_CohesiveElements',
                         sets=(rvePart.sets['set_AllElements'],
                               rvePart.sets['set_FiberElements'],
                               rvePart.sets['set_MatrixElements']),
                         operation=DIFFERENCE)

    elemType_coh = ElemType(elemCode=COH2D4, elemLibrary=STANDARD,
                            elemDeletion=ON, viscosity=0.0001, maxDegradation=0.99)

    if 'set_CohesiveElements' in rvePart.sets:
        coh_elements = rvePart.sets['set_CohesiveElements'].elements
        if coh_elements:
            rvePart.setElementType(regions=(coh_elements,), elemTypes=(elemType_coh,))
            rvePart.SectionAssignment(region=rvePart.sets['set_CohesiveElements'],
                                      sectionName='Section-Cohesive', offset=0.0)

    print("Step 5 Complete.")

    # ==================== 步骤 6: 周期性边界条件 ====================
    print("\nStep 6: Applying Periodic Boundary Conditions...")

    rootAssembly = a.rootAssembly
    rveInstance = rootAssembly.Instance(name='RVE-2D-1', part=rvePart, dependent=ON)

    createReferencePoints(a)

    dimensions = getRVEDimensions(a, 'RVE-2D-1')
    nodes_left, nodes_right, nodes_bottom, nodes_top = getBoundaryNodes(a, 'RVE-2D-1', dimensions)

    pairing_tolerance = globalSeedSize * pairingToleranceFactor

    if len(nodes_left) <= len(nodes_right):
        lr_pairs = pairBoundaryNodes(nodes_left, nodes_right, pairing_tolerance, 1)
    else:
        lr_pairs = [(p[1], p[0]) for p in pairBoundaryNodes(nodes_right, nodes_left, pairing_tolerance, 1)]

    if len(nodes_bottom) <= len(nodes_top):
        bt_pairs = pairBoundaryNodes(nodes_bottom, nodes_top, pairing_tolerance, 0)
    else:
        bt_pairs = [(p[1], p[0]) for p in pairBoundaryNodes(nodes_top, nodes_bottom, pairing_tolerance, 0)]

    print("  Left/Right: %d pairs from %d/%d nodes" %
          (len(lr_pairs), len(nodes_left), len(nodes_right)))
    print("  Bottom/Top: %d pairs from %d/%d nodes" %
          (len(bt_pairs), len(nodes_bottom), len(nodes_top)))

    applyPeriodicConstraints(a, 'RVE-2D-1', lr_pairs, 'Left-Right')
    applyPeriodicConstraints(a, 'RVE-2D-1', bt_pairs, 'Bottom-Top')

    print("Step 6 Complete.")

    print("\n" + "=" * 70)
    print("RVE MODEL GENERATION COMPLETED SUCCESSFULLY")
    print("=" * 70)


# =================================================================
#                 主程序入口
# =================================================================
if __name__ == '__main__':
    """主程序 - 设置参数并执行建模"""

    # ========== 核心几何参数 ==========
    # TARGET_VF: 目标纤维体积分数
    # - 取值范围: 0.0 ~ 0.7 (理论上可达0.785，但实际受最小间距限制)
    # - 说明: 定义纤维在RVE中所占的体积比例，根据此值自动计算纤维数量
    # - 影响: 值越大，复合材料强度越高，但生成难度增加
    TARGET_VF = 0.5

    # RVE_SIZE: RVE模型的尺寸 [宽度, 高度]，单位mm
    # - 说明: 代表性体积单元的物理尺寸
    # - 建议: 尺寸应至少为纤维直径的10倍以上，以保证统计代表性
    # - 影响: 尺寸越大，包含纤维数越多，计算成本越高，但统计性越好
    RVE_SIZE = [0.057, 0.057]

    # FIBER_RADIUS: 纤维半径，单位mm
    # - 说明: 单根纤维的半径，纤维直径 = 2 × 半径
    # - 影响: 决定了纤维的尺寸和数量（Vf相同时，半径越小纤维数越多）
    FIBER_RADIUS = 0.0035

    # ========== 约束参数 ==========
    # MIN_DIST_FACTOR: 最小间距因子
    # - 取值范围: 2.0 ~ 2.5 (推荐 2.0 ~ 2.1)
    # - 说明: 纤维中心之间的最小距离 = 因子 × 纤维半径
    #         例如：因子=2.05时，最小中心距=2.05R，纤维表面最小间隙=0.05R
    # - 影响:
    #   * 值越大: 纤维间隙越大，生成越容易，但可达到的最大Vf降低
    #   * 值越小: 可以达到更高的Vf，但生成难度大幅增加，可能失败
    # - 建议: 对于高Vf(>0.5)，使用2.0~2.05；对于低Vf(<0.4)，可用2.1~2.2
    MIN_DIST_FACTOR = 2.05

    # ========== 网格参数 ==========
    # GLOBAL_SEED_SIZE: 全局网格种子尺寸，单位mm
    # - 说明: 控制网格的精细程度，决定单元的平均边长
    # - 建议: 纤维半径的 1/7 ~ 1/10
    #         例如：半径0.0035mm时，种子尺寸0.0005mm (约1/7)
    # - 影响:
    #   * 值越小: 网格越精细，结果越准确，但计算成本显著增加
    #   * 值越大: 网格越粗糙，计算快但可能影响精度
    # - 注意: 纤维圆周至少需要20-30个单元以保证几何精度
    GLOBAL_SEED_SIZE = 0.0005

    # PAIRING_TOLERANCE_FACTOR: 边界节点配对容差因子
    # - 说明: 周期性边界条件中，左右/上下边界节点配对的容差
    #         实际容差 = 因子 × 网格尺寸
    # - 取值范围: 0.1 ~ 1.0 (推荐 0.5)
    # - 影响: 容差太小可能导致配对失败，太大可能误配对
    # - 建议: 一般使用默认值0.5即可，除非遇到配对问题
    PAIRING_TOLERANCE_FACTOR = 0.5

    # ========== 算法调优参数 ==========
    # RSA_SEEDING_RATIO: RSA播种比例
    # - 取值范围: 0.0 ~ 1.0
    # - 说明: 控制纤维分布模式的关键参数，决定有多少比例的纤维作为"锚点"
    # - 三种模式:
    #   * 高值 (0.8~1.0) - "均匀离散"模式:
    #     - 纤维均匀散布，避免团簇
    #     - 生成速度快（推荐用于生产环境）
    #     - 适合需要均匀分布的情况
    #   * 低值 (0.0~0.3) - "物理平衡"模式:
    #     - 模拟物理粒子排斥平衡
    #     - 会形成局部密集团簇（更接近真实材料）
    #     - 生成速度慢
    #   * 中等值 (0.4~0.7) - 混合模式:
    #     - 兼具两种模式的特点
    #     - 平衡速度和分布多样性
    # - 推荐: 0.9 (快速均匀分布) 或 0.1 (真实物理分布)
    RSA_SEEDING_RATIO = 0.9

    # ========== 纤维材料参数 (MPa) ==========
    # FIBER_E: 纤维弹性模量，单位MPa
    # - 说明: 纤维材料的杨氏模量，表征材料的刚度
    # - 典型值: 碳纤维 ~230000 MPa, 玻璃纤维 ~70000 MPa
    FIBER_E = 15000.0

    # FIBER_NU: 纤维泊松比
    # - 说明: 材料横向应变与轴向应变的比值
    # - 取值范围: 0.0 ~ 0.5 (大多数材料在 0.2 ~ 0.35)
    FIBER_NU = 0.3

    # ========== 基体材料参数 (MPa) ==========
    # MATRIX_E: 基体弹性模量，单位MPa
    # - 说明: 基体材料的杨氏模量
    # - 典型值: 环氧树脂 ~3000 MPa, 聚酯 ~3500 MPa
    MATRIX_E = 3170.0

    # MATRIX_NU: 基体泊松比
    # - 说明: 基体材料的泊松比
    # - 典型值: 聚合物基体 0.3 ~ 0.4
    MATRIX_NU = 0.35

    # MATRIX_FRICTION_ANGLE: Drucker-Prager摩擦角，单位度
    # - 说明: 描述材料的剪切强度随压应力变化的参数
    # - 取值范围: 0° ~ 45° (典型值 10° ~ 30°)
    # - 影响: 角度越大，材料受压时强度增加越明显
    MATRIX_FRICTION_ANGLE = 16.0

    # MATRIX_FLOW_STRESS_RATIO: Drucker-Prager流动应力比 K
    # - 说明: 屈服面在偏平面上的形状参数，K=1.0表示圆形屈服面
    # - 取值范围: 0.778 ~ 1.0
    # - 影响: 控制拉压屈服强度的差异
    MATRIX_FLOW_STRESS_RATIO = 1.0

    # MATRIX_DILATION_ANGLE: Drucker-Prager膨胀角，单位度
    # - 说明: 控制塑性流动时的体积变化
    # - 取值范围: 0° ~ 摩擦角
    # - 影响: 角度越大，塑性变形时体积膨胀越明显
    # - 建议: 通常取与摩擦角相同或略小的值
    MATRIX_DILATION_ANGLE = 16.0

    # MATRIX_HARDENING_YIELD: 硬化屈服应力，单位MPa
    # - 说明: 材料开始塑性硬化时的屈服应力
    # - 影响: 定义了材料从弹性到塑性的转变点
    MATRIX_HARDENING_YIELD = 106.4

    # MATRIX_HARDENING_PLASTIC_STRAIN: 硬化对应的塑性应变
    # - 说明: 与硬化屈服应力对应的塑性应变值
    # - 影响: 配合屈服应力定义硬化曲线
    # - 注意: 0.0表示初始屈服点
    MATRIX_HARDENING_PLASTIC_STRAIN = 0.0

    # MATRIX_DAMAGE_STRAIN: 韧性损伤起始应变
    # - 说明: 材料开始发生韧性损伤时的等效塑性应变
    # - 影响: 值越小，材料越早开始损伤
    # - 典型值: 0.01 ~ 0.1
    MATRIX_DAMAGE_STRAIN = 0.01

    # MATRIX_DAMAGE_STRESS_TRIAX: 应力三轴度
    # - 说明: 静水压力与等效应力的比值，η = σm / σeq
    # - 取值范围: -1/3 (纯剪) 到 +∞ (静水拉伸)
    # - 影响: 描述应力状态对损伤的影响
    # - 注意: 0.0表示不考虑三轴度影响
    MATRIX_DAMAGE_STRESS_TRIAX = 0.0

    # MATRIX_DAMAGE_STRAIN_RATE: 应变率，单位1/s
    # - 说明: 损伤发生时的应变率
    # - 影响: 描述应变率对损伤的影响
    # - 注意: 0.0表示准静态加载
    MATRIX_DAMAGE_STRAIN_RATE = 0.0

    # MATRIX_DAMAGE_DISPLACEMENT: 损伤演化特征位移，单位mm
    # - 说明: 从损伤起始到完全失效的特征位移
    # - 影响: 值越大，材料损伤演化越慢，韧性越好
    # - 典型值: 0.00001 ~ 0.0001 mm
    # - 注意: 应根据网格尺寸调整，建议为单元尺寸的0.1~1倍
    MATRIX_DAMAGE_DISPLACEMENT = 5e-05

    # ========== 界面材料参数 ==========
    # COHESIVE_K_NN: 界面法向刚度，单位N/mm³
    # - 说明: 控制界面法向（垂直于界面）的初始刚度
    # - 建议: 取10⁷ ~ 10⁹，通常为基体模量的10³~10⁵倍
    # - 影响: 刚度越大，界面初始响应越硬
    COHESIVE_K_NN = 1e8

    # COHESIVE_K_SS: 界面第一切向刚度，单位N/mm³
    # - 说明: 控制界面第一个切向的初始刚度
    # - 建议: 通常与法向刚度相同或略小
    COHESIVE_K_SS = 1e8

    # COHESIVE_K_TT: 界面第二切向刚度，单位N/mm³
    # - 说明: 控制界面第二个切向的初始刚度
    # - 建议: 2D模型中通常与K_SS相同
    COHESIVE_K_TT = 1e8

    # COHESIVE_T_N: 界面法向强度，单位MPa
    # - 说明: 界面在法向（拉伸）方向的最大承载能力
    # - 影响: 决定界面开裂的起始载荷
    # - 典型值: 20 ~ 80 MPa
    COHESIVE_T_N = 44.0

    # COHESIVE_T_S: 界面第一切向强度，单位MPa
    # - 说明: 界面在第一个切向（剪切）方向的最大承载能力
    # - 影响: 决定界面剪切失效的起始载荷
    # - 典型值: 40 ~ 100 MPa (通常大于法向强度)
    COHESIVE_T_S = 82.0

    # COHESIVE_T_T: 界面第二切向强度，单位MPa
    # - 说明: 界面在第二个切向方向的最大承载能力
    # - 建议: 2D模型中通常与T_S相同
    COHESIVE_T_T = 82.0

    # COHESIVE_GIC: I型断裂能（法向），单位N/mm
    # - 说明: 界面完全分离所需的能量（拉伸模式）
    # - 影响: 决定界面韧性和损伤演化速度
    # - 典型值: 0.0001 ~ 0.01 N/mm
    # - 注意: 值越大，界面越韧，失效过程越渐进
    COHESIVE_GIC = 0.001

    # COHESIVE_GIIC: II型断裂能（切向），单位N/mm
    # - 说明: 界面完全分离所需的能量（剪切模式）
    # - 影响: 决定剪切失效的韧性
    # - 典型值: 通常为GIC的1~3倍
    COHESIVE_GIIC = 0.002

    # COHESIVE_GIIIC: III型断裂能（切向），单位N/mm
    # - 说明: 界面完全分离所需的能量（撕裂模式）
    # - 建议: 2D模型中通常与GIIC相同
    COHESIVE_GIIIC = 0.002

    # COHESIVE_ETA: BK准则指数
    # - 说明: Benzeggagh-Kenane混合模式断裂准则的指数
    # - 取值范围: 1.0 ~ 3.0
    # - 影响: 控制混合模式下断裂能的插值方式
    # - 建议: 1.5 ~ 2.0 (材料试验拟合得到)
    COHESIVE_ETA = 1.5

    # COHESIVE_STAB_COEFF: 粘性稳定系数
    # - 说明: 用于数值稳定性的人工阻尼系数
    # - 取值范围: 0.00001 ~ 0.001
    # - 影响:
    #   * 值太小: 可能导致数值不稳定
    #   * 值太大: 可能影响结果准确性
    # - 建议: 0.0001 (需要时可微调)
    COHESIVE_STAB_COEFF = 0.0001

    # ========== 输出控制 ==========
    # EXPORT_COORDINATES: 是否导出纤维中心坐标到CSV文件
    # - True: 导出坐标，便于后续分析和验证
    # - False: 不导出
    EXPORT_COORDINATES = True

    # CSV_FILENAME: CSV文件名
    # - None: 自动生成文件名（格式：FiberCenters_Vf50_YYYYMMDD_HHMMSS.csv）
    # - 字符串: 使用指定的文件名
    CSV_FILENAME = None

    # ========== 模型命名 ==========
    targetModelName = '2D-RVE-Vf-%d' % (TARGET_VF * 100)
    tempModelName = targetModelName + '_TEMP_' + str(int(time.time()))

    print("\n" + "=" * 70)
    print("Starting 2D RVE Generation...")
    print("=" * 70)
    print("Target Model Name: %s" % targetModelName)

    if targetModelName in mdb.models:
        print("WARNING: Model exists and will be replaced.")

    print("Using temporary name: %s" % tempModelName)
    print("=" * 70)

    # ========== 执行建模 ==========
    createRVEModel(
        modelName=tempModelName,
        rveSize=RVE_SIZE,
        fiberRadius=FIBER_RADIUS,
        target_Vf=TARGET_VF,
        minDistanceFactor=MIN_DIST_FACTOR,
        globalSeedSize=GLOBAL_SEED_SIZE,
        pairingToleranceFactor=PAIRING_TOLERANCE_FACTOR,
        rsa_seeding_ratio=RSA_SEEDING_RATIO,
        export_coordinates=EXPORT_COORDINATES,
        csv_filename=CSV_FILENAME,
        fiber_E=FIBER_E,
        fiber_nu=FIBER_NU,
        matrix_E=MATRIX_E,
        matrix_nu=MATRIX_NU,
        matrix_friction_angle=MATRIX_FRICTION_ANGLE,
        matrix_flow_stress_ratio=MATRIX_FLOW_STRESS_RATIO,
        matrix_dilation_angle=MATRIX_DILATION_ANGLE,
        matrix_hardening_yield=MATRIX_HARDENING_YIELD,
        matrix_hardening_plastic_strain=MATRIX_HARDENING_PLASTIC_STRAIN,
        matrix_damage_strain=MATRIX_DAMAGE_STRAIN,
        matrix_damage_stress_triax=MATRIX_DAMAGE_STRESS_TRIAX,
        matrix_damage_strain_rate=MATRIX_DAMAGE_STRAIN_RATE,
        matrix_damage_displacement=MATRIX_DAMAGE_DISPLACEMENT,
        cohesive_K_nn=COHESIVE_K_NN,
        cohesive_K_ss=COHESIVE_K_SS,
        cohesive_K_tt=COHESIVE_K_TT,
        cohesive_t_n=COHESIVE_T_N,
        cohesive_t_s=COHESIVE_T_S,
        cohesive_t_t=COHESIVE_T_T,
        cohesive_GIC=COHESIVE_GIC,
        cohesive_GIIC=COHESIVE_GIIC,
        cohesive_GIIIC=COHESIVE_GIIIC,
        cohesive_eta=COHESIVE_ETA,
        cohesive_stab_coeff=COHESIVE_STAB_COEFF
    )

    # ========== 后处理 ==========
    if tempModelName in mdb.models:
        print("\n" + "=" * 70)
        print("Starting Model Cleanup and Rename...")
        print("=" * 70)

        models_to_delete = [m for m in mdb.models.keys()
                            if m != tempModelName]

        if models_to_delete:
            print("\nDeleting old models:")
            for m in models_to_delete:
                print("  - %s" % m)
            for modelKey in models_to_delete:
                del mdb.models[modelKey]
            print("Old models deleted.")
        else:
            print("\nNo old models to delete.")

        print("\nRenaming model:")
        print("  From: '%s'" % tempModelName)
        print("  To:   '%s'" % targetModelName)
        mdb.models.changeKey(fromName=tempModelName,
                             toName=targetModelName)

        print("\n" + "=" * 70)
        print("MODEL GENERATION COMPLETE")
        print("=" * 70)
    else:
        print("\n" + "!" * 70)
        print("ERROR: Model generation failed.")
        print("!" * 70 + "\n")