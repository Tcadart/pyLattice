# je laisse cela ici, pour le moment c'est la meme chose que points, mais si jamais a lavenir on aurait besoin d'un
# graphe la il faudrait separé entre points et node pour plus de clarté
class Node:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __eq__(self, other):
        return (
            isinstance(other, Node)
            and self.x == other.x
            and self.y == other.y
            and self.z == other.z
        )

    def __hash__(self):
        return hash((self.x, self.y, self.z))

    def __repr__(self):
        return f"Node({self.x}, {self.y}, {self.z})"
