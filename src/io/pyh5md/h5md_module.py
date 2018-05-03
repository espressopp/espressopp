import numpy as np
import h5py

CHUNK_0_MAX=32
CHUNK_1_MAX=128

class Element(object):
    def append(self, *args, **kwargs):
        raise NotImplementedError

class FixedElement(h5py.Dataset, Element):
    def __init__(self, loc, name, **kwargs):
        if name in loc:
            d = loc[name]
        else:
            d = loc.create_dataset(name, **kwargs)
        super(FixedElement, self).__init__(d._id)
        self.step = None
        self.step_offset = None
        self.time = None
        self.time_offset = None
    def append(self, v, step=None, time=None):
        pass
    def get_by_idx(self, idx):
        return self
    @property
    def value(self):
        return h5py.Dataset(self._id)
    @property
    def element_type(self):
        return 'FixedElement'
    def __repr__(self):
        return 'H5MD FixedElement'

class LinearElement(h5py.Group, Element):
    def __init__(self, loc, name, **kwargs):
        is_new = name not in loc
        g = loc.require_group(name)
        if is_new:
            step = kwargs.pop('step')
            step_offset = kwargs.pop('step_offset', None)
            time = kwargs.pop('time', None)
            time_offset = kwargs.pop('time_offset', None)
        super(LinearElement, self).__init__(g._id)
        if is_new:
            self.value = g.create_dataset('value', **kwargs)
            g.create_dataset('step', data=int(step))
            self.step = int(step)
            if step_offset is not None:
                g['step'].attrs['offset'] = int(step_offset)
                self.step_offset = int(step_offset)
            self.step_offset = None
            if time is not None:
                g.create_dataset('time', data=time)
                self.time = time
                if time_offset is not None:
                    g['time'].attrs['offset'] = time_offset
                self.time_offset = time_offset
        else:
            self.step = g['step'][()]
            if 'offset' in self['step'].attrs:
                self.step_offset = self['step'].attrs['offset']
            else:
                self.step_offset = None
            self.value = self['value']
            if 'time' in self:
                self.time = self['time'][()]
                if 'offset' in self['time'].attrs:
                    self.time_offset = self['time'].attrs['offset']
                else:
                    self.time_offset = None
            else:
                self.time = None
                self.time_offset = None
    @property
    def element_type(self):
        return 'LinearElement'
    def append(self, v, step=None, time=None, region=None, collective=False):
        self.value.resize(self.value.shape[0]+1, axis=0)
        if region is not None:
            if collective:
                with self.value.collective:
                    self.value[-1,region[0]:region[1],...] = v
            else:
                self.value[-1,region[0]:region[1],...] = v
        else:
            self.value[-1] = v
    def get_by_idx(self, idx):
        return self['value'][idx]
    def __repr__(self):
        return 'H5MD LinearElement'

class TimeElement(h5py.Group, Element):
    def __init__(self, loc, name, **kwargs):
        is_new = name not in loc
        g = loc.require_group(name)
        if is_new:
            step_from = kwargs.pop('step_from', None)
            if step_from is not None:
                g['step'] = step_from.step
                self.step = step_from.step
                self.own_step = False
            else:
                self.step = g.create_dataset('step', dtype=int, shape=(0,), maxshape=(None,))
                self.own_step = True
            time = kwargs.pop('time', None)
            if time is not None:
                if self.own_step:
                    if time:
                        self.time = g.create_dataset('time', dtype=float, shape=(0,), maxshape=(None,))
                    else:
                        raise ValueError("Time must be True or None for TimeElement")
                else:
                    g['time'] = self.time = step_from.time
            else:
                self.time = None
            self.value = g.create_dataset('value', **kwargs)
        else:
            self.step = g['step']
            if 'offset' in g['step'].attrs:
                self.step_offset = g['step'].attrs['offset']
            else:
                self.step_offset = None
            if 'time' in g:
                self.time = g['time']
                if 'offset' in g['time'].attrs:
                    self.time_offset = g['time'].attrs['offset']
                else:
                    self.time_offset = None
            else:
                self.time = None
                self.time_offset = None
            self.value = g['value']
        super(TimeElement, self).__init__(g._id)
    def append(self, v, step, time=None, region=None, collective=False):
        self.value.resize(self.value.shape[0]+1, axis=0)
        if region is not None:
            if collective:
                with self.value.collective:
                    self.value[-1,region[0]:region[1],...] = v
            else:
                self.value[-1,region[0]:region[1],...] = v
        else:
            self.value[-1] = v
        if self.own_step:
            self.step.resize(self.step.shape[0]+1, axis=0)
            self.step[-1] = step
            if self.time and len(self.time.shape)==1:
                self.time.resize(self.time.shape[0]+1, axis=0)
                self.time[-1] = time
    def get_by_idx(self, idx):
        return self['value'][idx]
    @property
    def element_type(self):
        return 'TimeElement'
    def __repr__(self):
        return 'H5MD TimeElement'

def default_chunks(shape):
    result = list(shape)
    if len(shape)==0:
        raise ValueError
    if shape==(0,):
        return (64,)
    if shape[0]==0:
        result[0] = CHUNK_0_MAX
        result[1] = min(CHUNK_1_MAX, shape[1])
        return tuple(result)
    else:
        result[:2] = CHUNK_0_MAX, CHUNK_1_MAX
        return tuple(np.minimum(result,shape))


def element(loc, name, **kwargs):
    if name in loc:
        tmp_element = loc[name]
        if isinstance(tmp_element,h5py.Group):
            assert 'value' in tmp_element
            assert 'step' in tmp_element
            if tmp_element['step'].shape==():
                return LinearElement(loc, name, **kwargs)
            else:
                return TimeElement(loc, name, **kwargs)
        elif isinstance(tmp_element, h5py.Dataset):
            return FixedElement(loc, name, **kwargs)
        else:
            return None
    store = kwargs.pop('store')
    if store=='fixed':
        return FixedElement(loc, name, **kwargs)
    if 'shape' in kwargs:
        assert 'data' not in kwargs
        kwargs['shape'] = (0,) + kwargs['shape']
        kwargs['chunks'] = default_chunks(kwargs['shape'])
        if 'maxshape' in kwargs:
            kwargs['maxshape'] = (None,) + kwargs['maxshape']
        else:
            kwargs['maxshape'] = (None,) + kwargs['shape'][1:]
    data = kwargs.pop('data', None)
    if data is not None:
        if type(data) is not np.ndarray:
            data = np.asarray(data)
        kwargs['shape'] = (0,) + data.shape
        if 'maxshape' in kwargs:
            kwargs['maxshape'] = (None,) + kwargs['maxshape']
        else:
            kwargs['maxshape'] = (None,) + data.shape
        kwargs['dtype'] = data.dtype
    if store=='linear':
        return LinearElement(loc, name, **kwargs)
    elif store=='time':
        return TimeElement(loc, name, **kwargs)
    else:
        raise ValueError
        
class File(h5py.File):
    def __init__(self, name, mode=None, *args, **kwargs):
        if mode=='w':
            author = kwargs.pop('author', 'N/A')
            author_email = kwargs.pop('author_email', None)
            creator = kwargs.pop('creator', 'N/A')
            creator_version = kwargs.pop('creator_version', 'N/A')
        super(File, self).__init__(name, mode, *args, **kwargs)
        if mode=='w':
            g = self.create_group('h5md')
            g.attrs['version'] = np.array([1,1])
            g = self.create_group('h5md/author')
            g.attrs['name'] = author
            if author_email is not None:
                g.attrs['email'] = author_email
            g = self.create_group('h5md/creator')
            g.attrs['name'] = creator
            g.attrs['version'] = creator_version
        else:
            if 'h5md' not in self:
                raise KeyError("h5md group not found in file")
            g = self['h5md']
            assert 'author' in g
            assert 'creator' in g

    def particles_group(self, name):
        return ParticlesGroup(self, name)

class ParticlesGroup(h5py.Group):
    """Represents a particles group within a H5MD file."""
    def __init__(self, parent, name):
        """Create a new ParticlesGroup object."""
        super(ParticlesGroup, self).__init__(parent.require_group('particles').require_group(name)._id)

    def create_box(self, dimension=None, boundary=None, **kwargs):
        self.box = self.require_group('box')
        if dimension is not None and boundary is not None:
            assert len(boundary)==dimension
            self.box.attrs['dimension'] = dimension
            self.box.attrs['boundary'] = np.string_(boundary)
        if len(kwargs)>0:
            self.box.edges = element(self.box, 'edges', **kwargs)
