import tkinter as tk
from tkinter import ttk
from tkinter import messagebox


def on_closing(win: tk.Tk) -> None:
    """
    обрабатываем закрытие программы
    запрашиваем подтверждение о закрытии
    и высылаем отчет.
    """
    if messagebox.askokcancel('Выход из приложения', 'Хотите выйти?'):
        win.destroy()


def create_window() -> tk.Tk:
    """
    Создаем главное окно
    :return: window
    """

    window = tk.Tk()
    window.title("RNBCoffee")
    w = window.winfo_screenwidth()
    h = window.winfo_screenheight()
    window.geometry("{}x{}".format(w-10, h-1))
    icon = tk.PhotoImage(file='logo.png')
    window.iconphoto(False, icon)
    return window


def run() -> None:

    win = create_window()
    win.protocol('WM_DELETE_WINDOW', lambda window=win: on_closing(window))
