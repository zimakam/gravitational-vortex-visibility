# gravitational-vortex-visibility
Gravity from vorticity via FFT-Poisson. ρ(x) = κ·Σ|ω(x)|², ∇²φ = 4πGρ. Novel model: mass emerges from vortex field. MIT.
# GravitationalVortex

> Gravity from vorticity via FFT-Poisson.
> ρ(x) = κ·Σ|ω(x)|², ∇²φ = 4πGρ — novel model where mass emerges from vortex field.

**Author:** Зиявутдинов Магомед Камалович (Zimaka)
**Email:** zimakam@gmail.com
**ORCID:** 0009-0005-9212-9921
**License:** MIT
**Version:** 1.0.0

---

## What is this

Vortex field is the source of gravity. Density is derived from vorticity:

ρ(x) = κ · Σ|ω_i(x)|²
∇²φ = 4πG·ρ
g = −∇φ

Unlike classical approach (mass → gravity), here gravity is a **consequence of vortex field**. Vorticity plays the role of mass.

## Uniqueness

### 1. Gravity from vorticity
Mass emerges from |ω|² — no need for separate mass distribution.

### 2. FFT-Poisson solver
Fast and exact solution via FFT.

### 3. Four vortex models
rankine, lamb_oseen, gaussian, fibonacci.

## Quick start

pip install numpy
python gravitational_vortex.py --selftest

## Author
- Зиявутдинов Магомед Камалович (Zimaka)
- ORCID: 0009-0005-9212-9921
- GitHub: @zimakam

## License
MIT