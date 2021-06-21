from .utils import namedtuple_from_dict
from .tl255 import TL255


F128 = namedtuple_from_dict(
    'F128',
    {
        # The F128 full Gaussian grid uses the same latitudes as TL255 (N128)
        'yvals': TL255.yvals,
    }
)
