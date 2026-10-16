import numpy as np
from math import sin, cos, radians
from datetime import datetime
from pathlib import Path
import matplotlib.pyplot as plt


# -----------------------------
# 1. 只用日期计算太阳赤纬（高精度近似）
# -----------------------------
def solar_declination(date):
    """基于年积日的高精度近似公式计算太阳赤纬（度）"""
    day = date.timetuple().tm_yday
    gamma = 2.0 * np.pi / 365.0 * (day - 1)
    delta_rad = (
        0.006918
        - 0.399912 * np.cos(gamma)
        + 0.070257 * np.sin(gamma)
        - 0.006758 * np.cos(2 * gamma)
        + 0.000907 * np.sin(2 * gamma)
        - 0.002697 * np.cos(3 * gamma)
        + 0.00148 * np.sin(3 * gamma)
    )
    return np.degrees(delta_rad)


# -----------------------------
# 2. 经纬度 -> 地心单位向量
# -----------------------------
def latlon_to_xyz(lat, lon):
    lat_r = radians(lat)
    lon_r = radians(lon)

    x = cos(lat_r) * cos(lon_r)
    y = cos(lat_r) * sin(lon_r)
    z = sin(lat_r)

    return np.array([x, y, z], dtype=float)


# -----------------------------
# 3. 求两个可能的太阳方向
# -----------------------------
def solve_sun_vectors(lat, lon, delta):
    """
    输入:
        lat, lon : 已知在晨昏线上的点
        delta    : 太阳赤纬（度）
    输出:
        两个可能的太阳方向向量（一个晨线，一个昏线）
    """
    r = latlon_to_xyz(lat, lon)

    sz = sin(radians(delta))
    a, b, rz = r[0], r[1], r[2]

    # a*sx + b*sy + rz*sz = 0
    # 且 sx^2 + sy^2 = 1 - sz^2
    c = -rz * sz
    rr = np.sqrt(a * a + b * b)

    if rr < 1e-12:
        raise ValueError("输入点过于接近两极，当前算法不适用")

    phi = np.arctan2(b, a)
    cos_term = c / rr

    if abs(cos_term) > 1:
        raise ValueError("该点在该日期不可能位于晨昏线上")

    alpha = np.arccos(np.clip(cos_term, -1.0, 1.0))
    rho = np.sqrt(max(0.0, 1.0 - sz * sz))

    theta1 = phi + alpha
    theta2 = phi - alpha

    sun_vectors = []
    for theta in [theta1, theta2]:
        sx = rho * np.cos(theta)
        sy = rho * np.sin(theta)
        s = np.array([sx, sy, sz], dtype=float)
        s = s / np.linalg.norm(s)
        sun_vectors.append(s)

    return sun_vectors


# -----------------------------
# 4. 自动判断晨线/昏线
# -----------------------------
def classify_terminator_at_point(lon, s):
    """
    在给定点处，自动判断该太阳方向对应晨线还是昏线

    原理:
      局地东向单位矢量 e = (-sin(lon), cos(lon), 0)
      若 s·e > 0，太阳在东方 -> 晨线
      若 s·e < 0，太阳在西方 -> 昏线
    """
    lon_r = radians(lon)
    east = np.array([-np.sin(lon_r), np.cos(lon_r), 0.0], dtype=float)

    east_component = float(np.dot(s, east))

    if east_component > 0:
        return "晨线"
    if east_component < 0:
        return "昏线"
    return "无法判定"


# -----------------------------
# 5. 根据太阳方向生成晨昏线（保证经过输入点）
# -----------------------------
def terminator_from_sun_vector_through_point(s, lat, lon, n=180):
    """
    晨昏线 = 平面(s·r=0) 与单位球的交线
    以输入点为参数起点，保证输出曲线经过输入点
    返回:
        lon, lat, xyz  (度, 度, 地心系单位向量)
    """
    s = s / np.linalg.norm(s)

    # 将输入点投影到晨昏平面，作为参数化起点。
    r0 = latlon_to_xyz(lat, lon)
    u = r0 - np.dot(r0, s) * s

    if np.linalg.norm(u) < 1e-10:
        # 极端退化时退回到通用基向量方式。
        base = np.array([0.0, 0.0, 1.0])
        if np.linalg.norm(np.cross(base, s)) < 1e-8:
            base = np.array([1.0, 0.0, 0.0])
        u = np.cross(s, base)

    u = u / np.linalg.norm(u)
    v = np.cross(s, u)
    v = v / np.linalg.norm(v)

    t = np.linspace(0, 2 * np.pi, n)
    circle = np.outer(np.cos(t), u) + np.outer(np.sin(t), v)

    x = circle[:, 0]
    y = circle[:, 1]
    z = circle[:, 2]

    term_lon = np.degrees(np.arctan2(y, x))
    term_lat = np.degrees(np.arcsin(np.clip(z, -1.0, 1.0)))

    # 强制把输入点写入曲线首尾，确保导出结果显式经过输入点。
    term_lon[0], term_lat[0] = lon, lat
    term_lon[-1], term_lat[-1] = lon, lat

    # 与强制写入的经纬度保持一致，首尾点地心坐标也显式写为输入点。
    r0 = latlon_to_xyz(lat, lon)
    circle[0] = r0
    circle[-1] = r0

    return term_lon, term_lat, circle


# -----------------------------
# 6. 主函数：日期 + 晨昏线上的点 -> 两条候选线，并自动标注晨/昏
# -----------------------------
def compute_terminators(date, lat, lon):
    delta = solar_declination(date)
    sun_vectors = solve_sun_vectors(lat, lon, delta)

    results = []
    for s in sun_vectors:
        kind = classify_terminator_at_point(lon, s)
        term_lon, term_lat, term_xyz = terminator_from_sun_vector_through_point(s, lat, lon)
        results.append(
            {
                "type": kind,
                "sun_vector": s,
                "terminator_lon": term_lon,
                "terminator_lat": term_lat,
                "terminator_xyz": term_xyz,
            }
        )

    # 为了输出稳定：晨线放前，昏线放后
    results.sort(key=lambda item: 0 if item["type"] == "晨线" else 1)
    return delta, results


# -----------------------------
# 7. 保存晨昏线经纬度到本地TXT
# -----------------------------

def save_terminators_to_txt(date, lat, lon, results, output_dir=None):
    if output_dir is None:
        output_dir = Path(__file__).resolve().parent
    else:
        output_dir = Path(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)
    date_str = date.strftime("%Y%m%d")
    out_path = output_dir / f"terminator_{date_str}_{lat:.2f}_{lon:.2f}.txt"

    with out_path.open("w", encoding="utf-8") as f:
        f.write(f"date={date.strftime('%Y-%m-%d')}\n")
        f.write(f"input_point=({lat:.6f}, {lon:.6f})\n")
        input_xyz = latlon_to_xyz(lat, lon)
        f.write(
            "input_point_xyz=(%.10f, %.10f, %.10f)\n"
            % (input_xyz[0], input_xyz[1], input_xyz[2])
        )
        f.write("columns: index, lon_deg, lat_deg, x, y, z\n\n")

        for item in results:
            f.write(f"[{item['type']}] points={len(item['terminator_lon'])}\n")
            for idx, (lo, la, xyz) in enumerate(
                zip(
                    item["terminator_lon"],
                    item["terminator_lat"],
                    item["terminator_xyz"],
                ),
                start=1,
            ):
                f.write(
                    f"{idx},{lo:.6f},{la:.6f},{xyz[0]:.10f},{xyz[1]:.10f},{xyz[2]:.10f}\n"
                )
            f.write("\n")

    return out_path


# -----------------------------
# 8. 示例
# -----------------------------
if __name__ == "__main__":
    date = datetime(2025, 3,20)
    lat = 19.5
    lon = 109.0

    delta, results = compute_terminators(date, lat, lon)

    print(f"日期: {date.strftime('%Y-%m-%d')}")
    print(f"输入点: ({lat}°, {lon}°)")
    print(f"太阳赤纬: {delta:.6f}°")
    print()

    for i, item in enumerate(results, 1):
        print(f"候选 {i}: {item['type']}")
        print("太阳方向向量 s =", np.round(item["sun_vector"], 8))
        print(f"经纬度点数量: {len(item['terminator_lon'])}")
        print(
            "首点(应接近输入点)=",
            f"({item['terminator_lat'][0]:.6f}°, {item['terminator_lon'][0]:.6f}°)",
        )
        print()

    txt_path = save_terminators_to_txt(date, lat, lon, results)
    print(f"晨昏线经纬度已保存到: {txt_path}")

    plt.figure(figsize=(10, 5))

    for item in results:
        plt.plot(item["terminator_lon"], item["terminator_lat"], label=item["type"])

    plt.scatter([lon], [lat], label="输入点")
    plt.xlim(-180, 180)
    plt.ylim(-90, 90)
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.title(f"Terminator on {date.strftime('%Y-%m-%d')}")
    plt.grid(True)
    plt.legend()
    plt.show()





