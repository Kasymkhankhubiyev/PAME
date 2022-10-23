import tkinter as tk
from tkinter import ttk


class MainPanel:
    """

    """
    def _draw_work_panel(self):
        pass

    def _draw_graph(self):
        pass

    def __init__(self, window: tk.Tk):
        self.win = window
        self.graph_canvas = tk.Canvas(self.win, height=700, width=700, bg='white')
        self.panel_canvas = tk.Canvas(self.win, height=700, width=300)

    def draw_window(self) -> None:

        self.graph_canvas.grid(row=0, column=0, padx=5, pady=5)
        self.panel_canvas.grid(row=0, column=1, pady=5, padx=5)

        # ttk.Separator(self.win, orient='horizontal').grid(row=1, column=0)  # place(x=720, y=725)

        font_size = 20

        tk.Label(self.panel_canvas, text='Эффективная масса', font=('Arial', font_size)).grid(row=0, column=0)
        self.effective_mass = tk.Spinbox(self.panel_canvas, width=7, font=font_size)
        self.effective_mass.grid(row=0, column=1)

        tk.Label(self.panel_canvas, text='Ширина ямы', font=('Arial', font_size)).grid(row=1, column=0)
        self.hole_width = tk.Spinbox(self.panel_canvas, width=7, font=font_size)
        self.hole_width.grid(row=1, column=1)

        tk.Label(self.panel_canvas, text='Глубина ямы', font=('Arial', font_size)).grid(row=2, column=0)
        self.hole_depth = tk.Spinbox(self.panel_canvas, width=7, font=font_size)
        self.hole_depth.grid(row=2, column=1)

        tk.Label(self.panel_canvas, text='Расстояние между ямами', font=('Arial', font_size)).grid(row=3, column=0)
        self.holes_dist = tk.Spinbox(self.panel_canvas, width=7, font=font_size)
        self.holes_dist.grid(row=3, column=1)

        tk.Label(self.panel_canvas, text='Уровень энергии', font=('Arial', font_size)).grid(row=4, column=0)
        self.energy_level = tk.Label(self.panel_canvas, text='', font=('Arial', font_size), width=7)
        self.energy_level.grid(row=5, column=1)




