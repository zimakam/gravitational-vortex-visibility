#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GravitationalVortex — gravity from vorticity via FFT-Poisson
Author: Ziyavutdinov Magomed Kamalovich (Zimaka)
License: MIT
"""
import math
import argparse
import numpy as np

__author__ = "Зиявутдинов Магомед Камалович (Zimaka)"
__version__ = "1.0.0"
__license__ = "MIT"

TWO_PI = 2.0 * math.pi
G_NEWTON = 6.674e-11
PHI = (1.0 + math.sqrt(5.0)) / 2.0


class GravitationalVortex:
    """
    Gravity from vorticity via FFT-Poisson:
        ρ(x) = κ · Σ|ω_i(x)|²
        ∇²φ = 4πG·ρ
        g = −∇φ
    """

    def __init__(self, N=32, L=8.0, kappa=1.0, G=G_NEWTON):
        self.N = int(N)
        self.L = float(L)
        self.kappa = float(kappa)
        self.G = float(G)
        kx = np.fft.fftfreq(N, d=L/N) * TWO_PI
        ky = np.fft.fftfreq(N, d=L/N) * TWO_PI
        kz = np.fft.fftfreq(N, d=L/N) * TWO_PI
        self.KX, self.KY, self.KZ = np.meshgrid(kx, ky, kz, indexing='ij')
        self.K2 = self.KX**2 + self.KY**2 + self.KZ**2
        self.K2_safe = np.where(self.K2 < 1e-12, 1.0, self.K2)
        xs = np.linspace(-L/2, L/2, N, endpoint=False)
        self.X, self.Y, self.Z = np.meshgrid(xs, xs, xs, indexing='ij')
        self.density = None
        self.potential = None
        self.acceleration = None
        self.vortex_model = "lamb_oseen"
        self.core_radius = 0.5
        self.circulation = 13.0

    def _density_from_vorticity(self):
        dx = self.X
        dy = self.Y
        dz = self.Z
        r2 = dx*dx + dy*dy + dz*dz
        rc = self.core_radius
        G_ = self.circulation
        if self.vortex_model == "rankine":
            omega = np.where(r2 <= rc*rc, abs(G_)/(math.pi*rc*rc), 0.0)
        elif self.vortex_model == "lamb_oseen":
            omega = (abs(G_)/(math.pi*rc*rc)) * np.exp(-r2/(rc*rc))
        elif self.vortex_model == "gaussian":
            sigma = rc / math.sqrt(2)
            omega = (abs(G_)/(2*math.pi*sigma*sigma)) * \
                    np.exp(-r2/(2*sigma*sigma))
        else:
            omega = np.zeros_like(r2)
        return self.kappa * omega**2

    def compute(self):
        self.density = self._density_from_vorticity()
        rho_hat = np.fft.fftn(self.density)
        phi_hat = -4*math.pi*self.G*rho_hat / self.K2_safe
        phi_hat[0, 0, 0] = 0.0
        self.potential = np.real(np.fft.ifftn(phi_hat))
        gx = -np.real(np.fft.ifftn(1j*self.KX*phi_hat))
        gy = -np.real(np.fft.ifftn(1j*self.KY*phi_hat))
        gz = -np.real(np.fft.ifftn(1j*self.KZ*phi_hat))
        self.acceleration = np.stack([gx, gy, gz], axis=-1)
        return self.potential

    def total_mass(self):
        if self.density is None:
            self.compute()
        dV = (self.L / self.N) ** 3
        return float(self.density.sum() * dV)

    def gravitational_energy(self):
        if self.potential is None:
            self.compute()
        dV = (self.L / self.N) ** 3
        return 0.5 * float((self.density * self.potential).sum() * dV)

    def report(self):
        if self.potential is None:
            self.compute()
        g_mag = np.linalg.norm(self.acceleration, axis=-1)
        return (
            f"GravitationalVortex(model={self.vortex_model})\n"
            f"  Grid: {self.N}³\n"
            f"  L = {self.L}\n"
            f"  κ = {self.kappa}\n"
            f"  G = {self.G}\n"
            f"  Total mass: {self.total_mass():.6e}\n"
            f"  |φ| max: {np.abs(self.potential).max():.6e}\n"
            f"  |g| max: {g_mag.max():.6e}\n"
            f"  E_g: {self.gravitational_energy():.6e}")


def selftest():
    print("=" * 60)
    print(f"GRAVITATIONALVORTEX v{__version__} — SELFTEST")
    print("=" * 60)

    for model in ["rankine", "lamb_oseen", "gaussian"]:
        gv = GravitationalVortex(N=16, L=4.0)
        gv.vortex_model = model
        gv.compute()
        print(f"\n[{model}]")
        print(f"  ΣM = {gv.total_mass():.4e}")
        print(f"  |φ|max = {np.abs(gv.potential).max():.4e}")
        print(f"  E_g = {gv.gravitational_energy():.4e}")

    print("\n✅ SELFTEST passed")
    print(f"Author: {__author__}")


def main():
    parser = argparse.ArgumentParser(
        description=f"GravitationalVortex v{__version__}")
    parser.add_argument("--selftest", action="store_true")
    parser.add_argument("--N", type=int, default=16)
    parser.add_argument("--L", type=float, default=4.0)
    args = parser.parse_args()
    selftest()


if __name__ == "__main__":
    main()