# -*- coding: utf-8 -*-
"""NeuroSpatialGradients (nsgrad)

TNBC 感觉神经驱动的空间免疫架构分型 —— 分析流水线包。

模块分工
--------
config     基因签名、空间尺度、QC 阈值、绘图风格 —— 全局唯一来源
io         Visium SpaceRanger 输出 -> AnnData(带微米坐标)
signature  基因 -> spot 级生物学读数
nerve      神经源点识别(三重过滤 + 空间连续性 + 阈值敏感性)
gradient   环形缓冲区 + 指数衰减拟合 y=a*exp(-b*d)+c + 离散降级方案
"""
__version__ = "0.1.0"

from . import config, io, signature, nerve, gradient  # noqa: F401
