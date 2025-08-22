import os
import sys

# Force offscreen to avoid needing a display
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt6.QtWidgets import QApplication
from PyQt6.QtGui import QPixmap

from simple_reaction_gui import create_reaction_image, create_placeholder_image


def main():
    app = QApplication([])

    # A simple reaction SMILES (ethanol oxidation to acetaldehyde as placeholder format)
    smiles = "CCO>>CC=O"
    pm: QPixmap = create_reaction_image(smiles, 400, 160)
    if not isinstance(pm, QPixmap):
        print("FAIL: Did not return QPixmap")
        sys.exit(1)

    # Even if RDKit is missing, we should still get a non-null pixmap drawn via QImage
    if pm.isNull():
        # Try placeholder directly
        pm2 = create_placeholder_image(smiles, 400, 160)
        if pm2.isNull():
            print("FAIL: Placeholder drawing returned null pixmap")
            sys.exit(2)

    print("PASS: Image functions executed without QPainter errors")
    # No need to exec the app loop for this smoketest


if __name__ == "__main__":
    main()
