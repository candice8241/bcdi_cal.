

import streamlit as st
st.set_page_config(page_title="FCC Bragg Calculator", layout="centered")
import numpy as np
st.title("📐 FCC Bragg & Coherence GUI")
st.header("📷 Detector + Sampling Parameters")
# -------------------
# Core Physics Functions
# -------------------
def bragg_theta(lambda_, d_spacing, n=1, deg=False):
    sin_theta = n * lambda_ / (2 * d_spacing)
    if sin_theta > 1:
        raise ValueError("无解：sin(theta) > 1，检查波长和晶面间距")
    theta = np.arcsin(sin_theta)
    return np.degrees(theta) if deg else theta

def d_fcc(a, h, k, l):
    return a / np.sqrt(h**2 + k**2 + l**2)

def compute_Bdet(delta, gamma):
    cosd, sind = np.cos(delta), np.sin(delta)
    cosg, sing = np.cos(gamma), np.sin(gamma)
    Bdet = np.array([
        [cosd, -sing*sind, cosg*sind],
        [0,     cosg,      sing],
        [-sind, -cosd*sing, cosd*cosg]
    ])
    return Bdet
def compute_Brecip(P, D, delta, gamma, delta_theta, lambda_):
    cosd, sind = np.cos(delta), np.sin(delta)
    cosg, sing = np.cos(gamma), np.sin(gamma)
    Brecip = np.array([
        [P/(lambda_*D)*cosd, -P/(lambda_*D)*sing*sind, delta_theta/lambda_*(1 - cosg*cosd)],
        [0, P/(lambda_*D)*cosg, 0],
        [-P/(lambda_*D)*sind, -P/(lambda_*D)*cosd*sing, delta_theta/lambda_*cosg*sind]
    ])
    return Brecip
def compute_Breal(Brecip, N1, N2, N3):
    D_inv = np.diag([1/N1, 1/N2, 1/N3])
    Breal = np.linalg.inv(Brecip.T) @ D_inv
    return Breal


# 输入区
E_keV = st.number_input("光子能量 (keV)", value=9.0, step=0.1)
a_angstrom = st.number_input("晶格常数 a (Å)", value=3.608, step=0.001)
h = st.number_input("h", value=1, step=1)
k = st.number_input("k", value=1, step=1)
l = st.number_input("l", value=1, step=1)
x_det = st.number_input("探测器像素尺寸 x_det (m)", value=55e-6, format="%.1e")
sigma = st.number_input("过采样率 σ", value=4.0, step=0.1)
sample_size = st.number_input("样品尺寸 (m)", value=1e-6, format="%.1e")
D = st.number_input("光斑尺寸 D (m)", value=100e-6, format="%.1e")
energy_resolution = st.number_input("能量分辨率 Δλ/λ", value=1e-4, format="%.1e")
space_resolution = st.number_input("空间分辨率 (m)", value=50e-9, format="%.1e")
N1 = st.number_input("N1", value=256, step=1)
N2 = st.number_input("N2", value=256, step=1)
N3 = st.number_input("N3", value=100, step=1)
gamma = st.number_input("Detector elevation γ (deg)", value=11.104)
delta = st.number_input("Detector azimuth  δ (deg)", value=29.607)

if st.button("计算"):
    try:
        # 转换 keV → λ
        lambda_ = 1.24e-9 / E_keV  

        # d_spacing
        d_spacing = d_fcc(a_angstrom * 1e-10, h, k, l)

        # Bragg angle
        delta_ = bragg_theta(lambda_, d_spacing, deg=True)

        # Calculations
        L_min = (x_det * sigma * sample_size) / lambda_
        delta_omega = lambda_ / (2 * np.sin(np.radians(delta_) / 2) * sample_size * sigma) * (180 / np.pi)
        L_T = (lambda_ / 2) * (L_min / D)
        L_L = lambda_ / (2 * energy_resolution)
        delta_q_speckle = (lambda_ * L_min) / sample_size
        q_max = 1 / (2 * space_resolution)
        delta_q = (2 * np.pi * x_det) / (L_min * lambda_)
        # 调用矩阵计算函数（确保单位角度转弧度）
        delta_rad = np.radians(delta)
        gamma_rad = np.radians(gamma)
        Bdet = compute_Bdet(delta_rad, gamma_rad)
        Brecip = compute_Brecip(x_det, D, delta_rad, gamma_rad, np.radians(delta_omega), lambda_)
        Breal = compute_Breal(Brecip, N1, N2, N3)* 1e9


        # 输出区
        st.subheader("📊 计算结果")
        st.write(f"λ = {lambda_:.3e} m")
        st.write(f"d(hkl) = {d_spacing:.3e} m")
        st.write(f"θ = {delta_:.3f} °")
        st.write(f"L_min = {L_min:.3f} m")
        st.write(f"Δω = {delta_omega:.4f} °")
        st.write(f"L_T (横向相干长度) = {L_T*1e6:.3f} µm")
        st.write(f"L_L (纵向相干长度) = {L_L*1e6:.3f} µm")
        st.write(f"Δq_speckle = {delta_q_speckle*1e6:.2f} µm")
        st.write(f"q_max = {q_max*1e-6:.2f} µm⁻¹")
        st.write(f"Δq = {delta_q*1e-9:.6f} nm⁻¹")
        st.text("B_det:")
        st.write(Bdet)
        st.text("B_recip: m⁻¹")
        st.write(Brecip)
        st.text("B_real: nm")
        st.write(Breal)                

    except Exception as e:
        st.error(f"错误: {e}")
