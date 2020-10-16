"""
FluorLabel PyMOL Plugin 

Label biomolecules with a fluorophore

(c) Fabio Steffen, University of Zurich, 2020
"""

dialog = None


def __init_plugin__(app=None):
    """
    Add FluorLabel plugin to the Plugins Menu
    """
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('FluorLabel', run_plugin_gui)


def run_plugin_gui():
    """
    Create the GUI Window
    """
    global dialog
    if dialog is None:
        from . import gui
        dialog = gui.App(_pymol_running=True)
    dialog.show()
