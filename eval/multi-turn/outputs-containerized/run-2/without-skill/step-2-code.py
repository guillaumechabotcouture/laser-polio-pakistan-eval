# Aggregate panel (bottom-right)
ax_agg = axes[1, 2]
ax_agg.plot(days, S_total / 1e3, label="S", color=COLORS["S"], lw=2)
ax_agg.plot(days, E_total / 1e3, label="E", color=COLORS["E"], lw=2)
ax_agg.plot(days, I_total / 1e3, label="I", color=COLORS["I"], lw=2)
ax_agg.plot(days, R_total / 1e3, label="R", color=COLORS["R"], lw=2)
ax_agg.axvline(peak_day, color="gray", ls="--", lw=1,
               label=f"I peak day {peak_day}")
ax_agg.set_title("All Patches Combined", fontsize=11)
ax_agg.set_xlabel("Day")
ax_agg.set_ylabel("Population (thousands)")
ax_agg.legend(loc="center right", fontsize=9)
ax_agg.grid(True, alpha=0.3)
ax_agg.set_xlim(0, NTICKS)

# Migration rate heatmap (top-right panel, previously unused)
ax_mig = axes[0, 2]
im = ax_mig.imshow(mig_rates * 100, cmap="YlOrRd", aspect="auto",
                   vmin=0, vmax=MAX_EXPORT_FRAC * 100)
ax_mig.set_title("Daily Migration Rates (%)", fontsize=11)
ax_mig.set_xticks(range(N_PATCHES))
ax_mig.set_yticks(range(N_PATCHES))
ax_mig.set_xticklabels([f"P{i}" for i in range(N_PATCHES)], fontsize=9)
ax_mig.set_yticklabels([f"P{i}" for i in range(N_PATCHES)], fontsize=9)
ax_mig.set_xlabel("Destination patch", fontsize=9)
ax_mig.set_ylabel("Origin patch", fontsize=9)
plt.colorbar(im, ax=ax_mig, label="%/day")
for i in range(N_PATCHES):
    for j in range(N_PATCHES):
        val = mig_rates[i, j] * 100
        txt_color = "white" if val > MAX_EXPORT_FRAC * 50 else "black"
        ax_mig.text(j, i, f"{val:.2f}", ha="center", va="center",
                    fontsize=8, color=txt_color)

plt.tight_layout()
plt.savefig("seir_4patch_gravity_results.png", dpi=150, bbox_inches="tight")
print("\nFigure saved: seir_4patch_gravity_results.png")
plt.show()