For a scalar SDE $dX(t) = f(X(t))\,dt + g(X(t))\,dW_t$, the Milstein update is
\begin{equation}
    X_{n+1} = X_n + f(X_n)\,\Delta t + g(X_n)\,\Delta W_n + \frac{1}{2}g(X_n)g'(X_n)\left((\Delta W_n)^2 - \Delta t\right).
\end{equation}
The only extra term compared with Euler--Maruyama is $\frac{1}{2}g(X_n)g'(X_n)\left((\Delta W_n)^2 - \Delta t\right)$, which only depends on $g(X)$.

\begin{itemize}
    \item \textbf{Deterministic ($g = 0$, Q4) and additive noise ($g = c$, Q2, Q3, Q5):} Here $g'(X) = 0$ regardless of whether $f(X)$ is $0$, $-\theta X$, or $X - X^3$. This results in the extra term vanishing, so the scheme reduces to Euler--Maruyama, offering no advantage.
    
    \item \textbf{Multiplicative noise ($f(X) = X - X^3$, $g(X) = cX$, Q6):} This case has an advantage, Since $g'(X) = c \neq 0$, the nonzero correction $\frac{1}{2}c^2 X_n\left((\Delta W_n)^2 - \Delta t\right)$ raises the strong order of convergence from $\gamma = \tfrac{1}{2}$ to $\gamma = 1$ (matching the weak order $\beta = 1$), giving a more accurate sample trajectory.
\end{itemize}
