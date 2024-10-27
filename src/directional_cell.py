from typing import Tuple


class DirectionalCell:
    """A class to represent a cell with 3 directions for backtracking"""
    def __init__(self, up: bool = False, diagonal: bool = False, left: bool = False):
        self.up = up
        self.diagonal = diagonal
        self.left = left
        self.num_directions = int(up) + int(diagonal) + int(left)

    def get_directions(self) -> Tuple[bool, bool, bool]:
        return self.up, self.diagonal, self.left
