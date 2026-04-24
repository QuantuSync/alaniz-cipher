from .field import FiniteField
from .graph import Graph
from .sheaf import Sheaf

# v3: extension field arithmetic for F_{p^d}
try:
    from .field_ext import ExtensionField
except ImportError:
    # galois dependency missing; v3 features unavailable
    pass
