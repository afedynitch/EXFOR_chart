import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import matplotlib.gridspec as gridspec

from .tools import (
    get_exfor_symb,
    elem_names,
    ExforDatabaseManager,
    NeucosDatabaseManager,
    periodictable,
    convert_to_full_tuple,
    data_dir,
    info
)

# Define drawing functions
xwidth = 0.9
ywidth = 0.9

colors = dict(
    endf="#000000",
    jendl="#ff9900",
    tendl="#808000",
    psb="#00ff00",
    neucosma="#666699",
    talys="#00ccff",
    peanut="#9933ff",
)


def gen_cds(from_list, db_name, xoffs=0.0, yoffs=0.0):
    return dict(
        db_name=[db_name for i in range(len(from_list))],
        nco_ids=[x[2] for x in from_list],
        A=[x[0] for x in from_list],
        N=[x[0] - x[1] for x in from_list],
        xpos=[x[0] - x[1] + xoffs for x in from_list],
        Z=[x[1] for x in from_list],
        ypos=[x[1] + yoffs for x in from_list],
        symx=[x[1] - 0.1 for x in from_list],
        elsymb=[x[3] for x in from_list],
    )


def draw_mpl_rect(
    axis, xcol, ycol, xwidth, ywidth, source, alpha, color, legend, picker=False
):
    nitems = len(list(source.items())[0][1])
    if nitems == 0:
        return (None, "_nolegend_")
    for i in range(len(list(source.items())[0][1])):
        axis.add_patch(
            patches.Rectangle(
                (source[xcol][i] - xwidth / 2.0, source[ycol][i] - ywidth / 2.0),
                xwidth,
                ywidth,
                color=color,
                alpha=alpha,
                picker=picker,
            )
        )

    proxy_artist = plt.Line2D(
        list(range(1)),
        list(range(1)),
        color="white",
        marker="s",
        markersize=14,
        alpha=alpha,
        markerfacecolor=color,
    )
    return (proxy_artist, legend)


def draw_mpl_cont(
    axis, xcol, ycol, xwidth, ywidth, source, alpha, color, legend, picker=False
):
    contour_width = 1.3
    nitems = len(list(source.items())[0][1])
    if nitems == 0:
        return (None, "_nolegend_")
    for i in range(len(list(source.items())[0][1])):
        axis.add_patch(
            patches.Rectangle(
                (source[xcol][i] - xwidth / 2.0, source[ycol][i] - ywidth / 2.0),
                xwidth,
                ywidth,
                edgecolor=color,
                alpha=alpha,
                facecolor="none",
                linewidth=contour_width,
                picker=picker,
            )
        )

    proxy_artist = plt.Line2D(
        list(range(1)),
        list(range(1)),
        color="white",
        marker="s",
        markersize=14,
        alpha=alpha,
        markerfacecolor="none",
        markeredgewidth=contour_width,
        markeredgecolor=color,
    )

    return (proxy_artist, legend)


def draw_mpl_circle(
    axis, xcol, ycol, radius, source, alpha, color, legend, picker=False
):
    nitems = len(list(source.items())[0][1])
    if nitems == 0:
        return (None, "_nolegend_")
    for i in range(len(list(source.items())[0][1])):
        axis.add_patch(
            patches.Circle(
                (source[xcol][i], source[ycol][i]),
                radius=radius,
                color=color,
                alpha=alpha,
                picker=picker,
            )
        )
    proxy_artist = plt.Line2D(
        list(range(1)),
        list(range(1)),
        color="white",
        marker="o",
        markersize=8,
        alpha=alpha,
        markerfacecolor=color,
    )
    return (proxy_artist, legend)


class PointBrowser(object):
    """
    Click on a point to select and highlight it -- the data that
    generated the point will be shown in the lower axes.  Use the 'n'
    and 'p' keys to browse through the next and previous points
    """

    def __init__(self, gs, fig, x4_dbm, nco_dbm, cds):
        self.lastind = 0
        self.x4_dbm, self.nco_dbm = x4_dbm, nco_dbm
        self.gs = gs
        self.fig = fig
        self.chart_ax = plt.subplot(gs[:2, :2])
        self.abs_ax = plt.subplot(gs[0, 2])
        self.n_ax = plt.subplot(gs[0, 3])
        self.p_ax = plt.subplot(gs[1, 2])
        self.n2_ax = plt.subplot(gs[1, 3])
        self.subpl_axes = [self.abs_ax, self.n_ax, self.p_ax, self.n2_ax]
        self.data = cds
        self.lastpick = (0, 0)
        self.lastshow = (0, 0)
        self.text = self.chart_ax.text(
            0.8,
            0.05,
            "",
            transform=self.chart_ax.transAxes,
            va="top",
            fontsize=18,
            weight="bold",
        )

        # self.selected, = ax.plot([xs[0]], [ys[0]], 'o', ms=12, alpha=0.4,
        #                          color='yellow', visible=False)

    def onpick(self, event):
        N = np.round(event.mouseevent.xdata).astype("int")
        Z = np.round(event.mouseevent.ydata).astype("int")
        nco_id = (N + Z) * 100 + Z
        if (N, Z) == self.lastpick:
            return
        self.lastpick = (N, Z)
        info(1, "{0}-{1} ({2})".format(elem_names[Z], N + Z, nco_id))
        # distances = np.hypot(x - xs[event.ind], y - ys[event.ind])
        # indmin = distances.argmin()
        # dataind = event.ind[indmin]

        self.lastind = nco_id
        info(1, "Updating.. from onpick")
        self.update()

    def onmove(self, event):
        if (
            event.xdata is None
            or event.ydata is None
            or event.xdata < -1.0
            or event.xdata > 30.0
            or event.ydata < -1
            or event.ydata > 30
        ):
            return
        N = np.round(event.xdata).astype("int")
        Z = np.round(event.ydata).astype("int")
        if (N, Z) == self.lastshow:
            return

        nco_id = (N + Z) * 100 + Z
        if N > 30 or Z > 30 or N < 0 or Z < 0:
            return
        self.text.set_text("%s" % get_exfor_symb(nco_id))
        self.fig.canvas.draw()
        self.lastshow = (N, Z)

    def draw_model_comparison(self, nco_id):
        for reac, ax in zip(["ABS", "N", "P", "2N"], self.subpl_axes):
            x4_dbm.set_reaction(reac)
            select = x4_dbm.get_entries(nco_id)
            for k, e in select.items():
                e.plot(fmt="s", label="short", axes=ax)

        for db, label, c in [
            (self.nco_dbm.neucos, "TALYS (NCO)", colors["neucosma"]),
            (self.nco_dbm.talys, "TALYS def.", colors["talys"]),
            (self.nco_dbm.psb, "PSB", colors["psb"]),
            (self.nco_dbm.endf, "ENDL-B-VII.1", colors["endf"]),
            (self.nco_dbm.jendl, "JENDL-PD-2004", colors["jendl"]),
            (self.nco_dbm.peanut, "PEANUT", colors["peanut"]),
            # (self.nco_dbm.tendl, 'TENDL-2015', colors['tendl']),
        ]:
            db.plot(
                (0, nco_id), label=label, lw=2, color=c, alpha=0.9, axes=self.abs_ax
            )
            db.plot(
                (nco_id, nco_id - 100),
                label=label,
                lw=2,
                color=c,
                alpha=0.9,
                axes=self.n_ax,
            )
            db.plot(
                (nco_id, nco_id - 101),
                label=label,
                lw=2,
                color=c,
                alpha=0.9,
                axes=self.p_ax,
            )
            db.plot(
                (nco_id, nco_id - 200),
                label=label,
                lw=2,
                color=c,
                alpha=0.9,
                axes=self.n2_ax,
            )

        titles = [
            "{0}(g,abs)".format(get_exfor_symb(nco_id)),
            "{0}(g,n){1}".format(get_exfor_symb(nco_id), get_exfor_symb(nco_id - 100)),
            "{0}(g,p){1}".format(get_exfor_symb(nco_id), get_exfor_symb(nco_id - 101)),
            "{0}(g,2n){1}".format(get_exfor_symb(nco_id), get_exfor_symb(nco_id - 200)),
        ]
        for t, ax in zip(titles, self.subpl_axes):
            ax.text(0.03, 0.97, t, transform=ax.transAxes, va="top")
            ax.set_xlabel("Photon energy in MeV")
            ax.set_ylabel("Cross-section in mb")
            # ax.set_title(t)
            ax.set_xlim(0, 80)
            ax.legend(loc="upper right", frameon=False, fontsize=12, numpoints=1)

    def update(self):
        if self.lastind is None:
            return

        dataind = self.lastind
        for ax in self.subpl_axes:
            ax.cla()
        # self.ax2.text(0.05, 0.9, 'Details for {0}'.format(
        #               get_exfor_symb(dataind)),
        #               transform=ax2.transAxes, va='top')
        self.draw_model_comparison(dataind)
        # self.selected.set_visible(True)
        # self.selected.set_data(xs[dataind], ys[dataind])

        self.text.set_text("%s" % get_exfor_symb(dataind))
        self.fig.canvas.draw()


x4_dbm = ExforDatabaseManager(debug=0)
nco_dbm = NeucosDatabaseManager()


# Retrieve data from EXFOR database
x4_dbm.set_reaction("ABS")
# set filters accordingly
# x4_dbm.set_filter('EVAL',False)
select = x4_dbm.get_entries(-1)


# Collect isotopes for which TALYS cross-sections are calculated
talys_requested = []
for el in periodictable.elements:
    if el.number > 26:
        break
    for iso in el.isotopes:
        exfor_symb = "{0}-{1}".format(el.symbol.upper(), iso)
        nco_num = 100 * iso + el.number
        talys_requested.append((iso, el.number, nco_num, exfor_symb))


# Collect isotopes from the reduced scheme for activation threshold 0.01
# (see note.pdf Fig. 20 center)
with open(data_dir / "neucosma_reduced_0_01.dat") as f:
    fline = f.readline().strip()
el_ids = [int(el.strip()) for el in fline.split(" ") if int(el.strip()) != 0]
nco_reduced = []
for elid in sorted(el_ids):
    Z = elid % 100
    A = elid / 100
    exfor_symb = "{0}-{1}".format(elem_names[Z].upper(), A)
    nco_reduced.append((A, Z, elid, exfor_symb))


seletion_boxes = []
for A in range(1, 30):
    for Z in range(1, 30):
        seletion_boxes.append(A * 100 + Z)

# Create "ColumnDataSource" tuples for convinient plotting (like bokeh)
selection_cds = gen_cds(convert_to_full_tuple(seletion_boxes), "all")
# nco_red_cds = gen_cds(nco_reduced, 'NCO, reduced')
nco_red_cds = gen_cds(convert_to_full_tuple(nco_dbm.neucos_abs_elems), "NCO")
exfor_cds = gen_cds(
    convert_to_full_tuple(x4_dbm.get_nco_list_of_selection()), "EXFOR, ABS"
)
exfor_all_cds = gen_cds(
    convert_to_full_tuple(list(set(np.loadtxt(data_dir / "exfor_all.txt")))),
    "EXFOR, *",
)

# Collect tables from WW calculations
try:
    ww_important = pickle.load(open(data_dir / "ww_important_isotopes.ppd", "rb"))
    nco_inj_cds = gen_cds(convert_to_full_tuple(ww_important[0]), "NCO, injected")
    nco1_cds = gen_cds(convert_to_full_tuple(ww_important[1]), "NCO, Level 1")
    nco2_cds = gen_cds(convert_to_full_tuple(ww_important[2]), "NCO, Level 2")
except FileNotFoundError:
    import warnings
    warnings.warn("NCO isotope collection not found, skipping.")
    nco_inj_cds = gen_cds([], "NCO, injected")
    nco1_cds = gen_cds([], "NCO, Level 1")
    nco2_cds = gen_cds([], "NCO, Level 2")

endf_cds = gen_cds(
    convert_to_full_tuple(nco_dbm.endf_abs_elems), "ENDF-B-VII.1", xoffs=-0.3, yoffs=0.3
)
jendl_cds = gen_cds(
    convert_to_full_tuple(nco_dbm.jendl_abs_elems),
    "JENDL-PD-2004",
    xoffs=+0.3,
    yoffs=0.3,
)
# tendl_cds = gen_cds(convert_to_full_tuple(nco_dbm.tendl_abs_elems),
#                     'TENDL-2015', xoffs=+0.3, yoffs=.3)
psb_cds = gen_cds(
    convert_to_full_tuple(nco_dbm.psb_abs_elems), "PSB", xoffs=-0.3, yoffs=-0.3
)
peanut_cds = gen_cds(
    convert_to_full_tuple(nco_dbm.peanut_abs_elems), "PEANUT", xoffs=0.3, yoffs=-0.3
)
peanut_cont_cds = gen_cds(
    convert_to_full_tuple(
        [el for el in nco_dbm.peanut_abs_elems if el not in nco_dbm.neucos_abs_elems]
    ),
    "_nolegend_",
)


def draw_iso_chart(chart_ax, logo=True):
    # chart_ax.set_title("Photo-nuclear cross-sections")

    leg_list = []
    leg_list.append(
        draw_mpl_rect(
            chart_ax,
            "xpos",
            "ypos",
            xwidth,
            ywidth,
            source=selection_cds,
            alpha=0.0,
            color="w",
            legend="_nolegend_",
            picker=True,
        )
    )
    leg_list.append(
        draw_mpl_rect(
            chart_ax,
            "xpos",
            "ypos",
            xwidth,
            ywidth,
            source=nco_red_cds,
            alpha=1.0,
            color="#bfbfbf",
            legend="TALYS 1.8 + CRPropa 2",
        )
    )
    leg_list.append(
        draw_mpl_rect(
            chart_ax,
            "xpos",
            "ypos",
            xwidth,
            ywidth,
            source=nco2_cds,
            alpha=0.6,
            color="#99ddff",
            legend="Astro, priority 2",
        )
    )
    leg_list.append(
        draw_mpl_rect(
            chart_ax,
            "xpos",
            "ypos",
            xwidth,
            ywidth,
            source=nco1_cds,
            alpha=0.6,
            color="#3333ff",
            legend="Astro, priority 1",
        )
    )
    leg_list.append(
        draw_mpl_rect(
            chart_ax,
            "xpos",
            "ypos",
            xwidth,
            ywidth,
            source=exfor_all_cds,
            alpha=0.9,
            color="#ffff66",
            legend=r"EXFOR, any",
        )
    )
    leg_list.append(
        draw_mpl_rect(
            chart_ax,
            "xpos",
            "ypos",
            xwidth,
            ywidth,
            source=exfor_cds,
            alpha=0.9,
            color="#ff0066",
            legend=r"EXFOR, $\sigma_{\rm abs}$",
        )
    )

    leg_list.append(
        draw_mpl_circle(
            chart_ax,
            "xpos",
            "ypos",
            0.13,
            source=endf_cds,
            alpha=0.9,
            color=colors["endf"],
            legend="ENDF-B-VII.1",
        )
    )
    leg_list.append(
        draw_mpl_circle(
            chart_ax,
            "xpos",
            "ypos",
            0.13,
            source=jendl_cds,
            alpha=0.9,
            color=colors["jendl"],
            legend="JENDL-PD-2004",
        )
    )
    # leg_list.append(draw_mpl_circle(chart_ax, "xpos", "ypos", 0.13, source=tendl_cds,
    #                 alpha=0.9, color=colors['tendl'], legend='TENDL-2015'))
    leg_list.append(
        draw_mpl_circle(
            chart_ax,
            "xpos",
            "ypos",
            0.13,
            source=psb_cds,
            alpha=0.8,
            color=colors["psb"],
            legend="PSB",
        )
    )
    leg_list.append(
        draw_mpl_circle(
            chart_ax,
            "xpos",
            "ypos",
            0.13,
            source=peanut_cds,
            alpha=0.8,
            color=colors["peanut"],
            legend="PEANUT",
        )
    )
    leg_list.append(
        draw_mpl_cont(
            chart_ax,
            "xpos",
            "ypos",
            xwidth + 0.02,
            ywidth + 0.02,
            source=peanut_cont_cds,
            alpha=0.3,
            color=colors["peanut"],
            legend="_nolegend_",
        )
    )
    leg_list.append(
        draw_mpl_cont(
            chart_ax,
            "xpos",
            "ypos",
            xwidth + 0.02,
            ywidth + 0.02,
            source=nco_inj_cds,
            alpha=1.0,
            color="#000000",
            legend="Astro, injected",
        )
    )
    leg_list = [leg for leg in leg_list if leg[1] != "_nolegend_"]
    chart_ax.set_xlim(-1, 31)
    chart_ax.set_ylim(-1, 27)
    chart_ax.set_xlabel
    chart_ax.set_xlabel("N")
    chart_ax.set_ylabel("Z")
    chart_ax.xaxis.label.set_fontsize(14)
    chart_ax.yaxis.label.set_fontsize(14)
    for label in chart_ax.get_xticklabels() + chart_ax.get_yticklabels():
        label.set_fontsize(14)
    chart_ax.grid()
    chart_ax.set_aspect("equal", "datalim")
    chart_ax.legend(*zip(*leg_list), loc="upper left", numpoints=1, frameon=False)
    if logo:
        import matplotlib.image as image

        nco_logo = image.imread("../files/neucos-logo-noborder.png")
        chart_ax.imshow(
            nco_logo, aspect="auto", extent=(25.2, 29.8, 0.5, 3.3), alpha=0.7, zorder=-1
        )


def start_interactive(show_logo):
    fig = plt.figure(figsize=(17, 7.5))
    gs = gridspec.GridSpec(2, 4)
    chart_ax = plt.subplot(gs[:2, :2])
    draw_iso_chart(chart_ax, show_logo)
    browser = PointBrowser(gs, fig, x4_dbm, nco_dbm, nco_red_cds)

    fig.canvas.mpl_connect("pick_event", browser.onpick)
    fig.canvas.mpl_connect("motion_notify_event", browser.onmove)
    plt.tight_layout()
    plt.show()


def start_static(show_logo, save):
    plt.figure(figsize=(8.5, 7.5))
    chart_ax = plt.subplot(111)
    draw_iso_chart(chart_ax, show_logo)

    plt.tight_layout()
    if save:
        plt.savefig("exfor_chart.pdf", dpi=300)
        plt.savefig("exfor_chart.png", dpi=300)
    plt.show()
