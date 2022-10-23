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

        # self.effective_mass = tk.Spinbox(self.win, width=7)
        # self.effective_mass.place(x=760, y=50)




