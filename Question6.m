\begin{figure}[htbp]
    \centering
    \includegraphics[width=\textwidth]{images/q6_multiplicative_noise.jpg}
    \caption{\textbf{Multiplicative noise dynamics ($g(X) = cX$).} \textit{Left:} Sample trajectories from $X(0) = 0.5$ for $c \in \{0.4, 1.1, 1.8\}$ compared with the deterministic equilibrium $X^* = 1$. \textit{Right:} Simulated stationary density ($M = 10^4$ paths, $c = 0.6$) overlaid with the analytical Fokker--Planck curve $p_s(x)$.}
    \label{fig:q6_multiplicative}
\end{figure}

\subsubsection*{How these solutions compare with previous types of noise}
\begin{itemize}
    \item \textbf{Elimination of Jumping Between Stable States:} Unlike additive noise ($g(X) = c$), where constant fluctuations push trajectories back and forth across $X = 0$ between the potential wells, multiplicative noise vanishes at the origin ($g(0) = 0$). Consequently, $X = 0$ acts as an absorbing barrier, ensuring that trajectories starting at $X(0) > 0$ remain strictly positive for all $t > 0$ and preventing switching between the $\pm 1$ equilibria.
    \item \textbf{Noise-Induced Stabilisation vs.\ Broadening:} Under additive noise ($g(X) = c$), increasing $c$ spreads the trajectories around $\pm 1$ and makes them jump more often between the two stable states at $X^* = +1$ and $X^* = -1$. Under multiplicative noise ($g(X) = cX$), increasing $c$ above a critical threshold ($|c| > \sqrt{2}$) results in all paths collapsing onto $X = 0$.

\end{itemize}

\subsubsection*{Key qualitative differences in the solutions}
\begin{itemize}
    \item \textbf{Panel 1 (Stochastic Stabilisation Threshold at $|c| = \sqrt{2}$):} Applying It\^o's lemma to $Y(t) = \ln X(t)$ for $X(t) > 0$ yields
    \begin{equation}
        d(\ln X) = \frac{1}{X}\,dX - \frac{1}{2X^2}(dX)^2 = \left(1 - X^2 - \frac{c^2}{2}\right)dt + c\,dW_t.
    \end{equation}
    When $X \approx 0$, there is an effective exponential growth rate of $\lambda = 1 - \frac{c^2}{2}$. For $|c| < \sqrt{2}$, $\lambda > 0$, so $X = 0$ remains unstable: the paths fluctuate around $X^* = 1$ under weak noise and drop to near zero when $1 \le c < \sqrt{2}$. Once $|c| > \sqrt{2}$, the negative noise $-\frac{c^2}{2}$ outweighs the $+1$ growth, pulling all paths towards $X = 0$, stabilising the origin.
    
    \item \textbf{Panel 2 (Skewed Stationary Density and Mode Shift):} For $0 < |c| < 1$ and $X(0) > 0$, solving the stationary Fokker--Planck equation gives a one-sided, right-skewed density. Differentiating $p_s(x)$ shows that the mode sits at  $x_{\mathrm{mode}} = \sqrt{1 - c^2}$, resulting in $x_{\mathrm{mode}} = 0.8$ for $c = 0.6$. Since the noise term grows away from the origin, the histogram peak shifts inward compared with $X^* = 1$.
    \end{itemize}



%figure

rng(42); clf;
dt = 0.005; T_traj = 30; T_stat = 50;
N_traj = round(T_traj / dt);
N_stat = round(T_stat / dt);
t_grid = linspace(0, T_traj, N_traj + 1);

figure('Position', [100, 100, 1000, 400]);

subplot(1, 2, 1); hold on;
c_vals = [0.4, 1.1, 1.8];
colors = {'#0072BD', '#D95319', '#7E2F8E'};
regimes = {'weak', 'intermittent', 'stabilised'};

for k = 1:length(c_vals)
    c = c_vals(k);
    X = zeros(1, N_traj + 1);
    X(1) = 0.5;
    dW = sqrt(dt) * randn(1, N_traj);
    for n = 1:N_traj
        X(n+1) = X(n) + (X(n) - X(n)^3)*dt + c*X(n)*dW(n) ...
            + 0.5*(c^2)*X(n)*(dW(n)^2 - dt);
    end
    plot(t_grid, X, 'Color', colors{k}, 'LineWidth', 1.2, ...
        'DisplayName', sprintf('c = %.1f (%s)', c, regimes{k}));
end
yline(1, '--', 'Color', '#EDB120', 'LineWidth', 1.3, 'HandleVisibility', 'off');
text(30.4, 1, 'X^* = 1', 'Color', '#EDB120', 'FontWeight', 'bold', 'FontSize', 10, ...
    'VerticalAlignment', 'middle', 'Clipping', 'off');
xlabel('Time'); ylabel('X(t)');
legend('Location', 'northeast'); grid on; box on;

subplot(1, 2, 2); hold on;
c_stat = 0.6; M = 10000;
X_ens = 0.5 * ones(M, 1);

for n = 1:N_stat
    dW = sqrt(dt) * randn(M, 1);
    X_ens = X_ens + (X_ens - X_ens.^3)*dt + c_stat*X_ens.*dW ...
        + 0.5*(c_stat^2)*X_ens.*(dW.^2 - dt);
end

histogram(X_ens, 60, 'Normalization', 'pdf', ...
    'FaceColor', [0.3 0.6 0.8], 'EdgeColor', 'none', ...
    'DisplayName', 'Empirical PDF (T = 50)');

x_pdf = linspace(0.01, 2.2, 400);
unnorm_pdf = x_pdf.^(2/c_stat^2 - 2) .* exp(-x_pdf.^2 / c_stat^2);
norm_pdf = unnorm_pdf / trapz(x_pdf, unnorm_pdf);
plot(x_pdf, norm_pdf, 'r-', 'LineWidth', 2, 'DisplayName', 'Analytical p_s(x)');
xline(sqrt(1 - c_stat^2), '--', 'Color', '#EDB120', 'LineWidth', 1.3, 'HandleVisibility', 'off');
text(sqrt(1 - c_stat^2), 1.53, 'Mode = 0.8', 'Color', '#EDB120', 'FontWeight', 'bold', 'FontSize', 10, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Clipping', 'off');

xlabel('x'); ylabel('Probability Density');
legend('Location', 'northeast'); grid on; box on;

exportgraphics(gcf, fullfile(matlabdrive, 'q6_multiplicative_noise.jpg'), 'Resolution', 300);
