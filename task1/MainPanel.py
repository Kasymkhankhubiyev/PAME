import tkinter as tk
from tkinter import ttk


class MainPanel:
    """

    """

    def __init__(self, window: tk.Tk):
        self.win = window
        self.canvas = tk.Canvas(self.win, height=700, width=700, bg='white')

    def draw_window(self) -> None:

        self.canvas.place(x=10, y=5)



