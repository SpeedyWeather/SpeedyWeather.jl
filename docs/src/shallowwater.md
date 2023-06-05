# Shallow water model

```math
\begin{aligned}
\frac{\partial \zeta}{\partial t} + \nabla \cdot (\mathbf{u}(\zeta + f)) &= (-1)^{n+1}\nu\nabla^{2n}\zeta, \\
\frac{\partial \mathcal{D}}{\partial t} - \nabla \times (\mathbf{u}(\zeta + f)) &= -\nabla^2(\tfrac{1}{2}(u^2 + v^2) + g\eta) + (-1)^{n+1}\nu\nabla^{2n}\mathcal{D}, \\
\frac{\partial \eta}{\partial t} + \nabla \cdot (\mathbf{u}h) &= 0,
\end{aligned}
```

where ``\zeta = \hat{\mathbf{z}} \cdot (\nabla \times \mathbf{u})`` is the relative vorticity,
``\mathcal{D} = \nabla \cdot \mathbf{u}`` the divergence, and ``\eta`` the deviation from the
fluid's rest height.

**Note**: more to come...

### Scaled shallow water equations

Similar to the scaled barotropic vorticity equations, the scaled shallow water equations scale the vorticity and the divergence equation with ``R^2``, but the continuity equation with ``R``

```math
\begin{aligned}
\frac{\partial \tilde{\zeta}}{\partial \tilde{t}} + \tilde{\nabla} \cdot (\mathbf{u}(\tilde{\zeta} + \tilde{f})) &=
\tilde{\nu}\tilde{\nabla}^{2n}\tilde{\zeta} \\
\frac{\partial \tilde{\mathcal{D}}}{\partial \tilde{t}} - \tilde{\nabla} \times (\mathbf{u}(\tilde{\zeta} + \tilde{f})) &=
-\tilde{\nabla}^2\left(\tfrac{1}{2}(u^2 + v^2) + g\eta \right) + \tilde{\nu}\tilde{\nabla}^{2n}\tilde{\mathcal{D}} \\
\frac{\partial \eta}{\partial \tilde{t}} + \tilde{\nabla} \cdot (\mathbf{u}h) &= 0.
\end{aligned}
```

## References

[^1]: Geophysical Fluid Dynamics Laboratory, [The Shallow Water Equations](https://www.gfdl.noaa.gov/wp-content/uploads/files/user_files/pjp/shallow.pdf).
[^2]: Geophysical Fluid Dynamics Laboratory, [The Spectral Dynamical Core](https://www.gfdl.noaa.gov/wp-content/uploads/files/user_files/pjp/spectral_core.pdf)