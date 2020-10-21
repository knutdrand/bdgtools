import numpy as np
import logging
from collections import Counter
from itertools import chain
from more_itertools import pairwise
from .regions import Regions

def broadcast(values, offsets):
    broadcasted = np.zeros(offsets[-1], dtype=values.dtype)
    broadcasted[offsets[1:-1]] = np.diff(values)
    broadcasted[0] = values[0]
    return np.cumsum(broadcasted)

class BedGraph:
    def __init__(self, indices, values, size=None, strict=True):
        self._indices = np.asanyarray(indices)
        assert np.issubdtype(self._indices.dtype, np.integer), self._indices
        if size is not None:
            assert np.all(self._indices<size), (self._indices, size)

        self._values = np.asanyarray(values)
        if size is not None:
            self._size = int(size)
        else:
            self._size = None

    def __iter__(self):
        pairs = zip(self._indices, chain(self._indices[1:], [self._size]))
        return zip(pairs, self._values)

    def _getitem(self, index):
        idx = np.searchsorted(self._indices, index, side="right")-1
        if idx >= 0:
            return self._values[idx]
        else:
            return 0

    def sum(self):
        if self._size is not None:
            sizes = np.diff(self._indices, append=self._size)
            return np.sum(sizes*self._values)
        assert self._values[-1] == 0, self._values[-1]
        return np.diff(self._indices*self._values[:-1])

    def mean(self):
        assert self._indices[0]==0
        return self.sum()/self._size

    def hist(self):
        h = np.zeros(int(np.max(self._values))+1)
        for (start, end), value in self:
            h[int(value)]+=(end-start)
        return h

    def reverse(self):
        assert self._size is not None
        indices = self._size-self._indices[::-1]
        assert indices[0] != 0
        indices = np.insert(indices[:-1], 0, 0)
        values = self._values[::-1]
        return self.__class__(indices, values, self._size)

    def _get_slice_indexes(self, start_idxs, end_idxs, directions, offsets):
        all_directions = broadcast(directions, offsets)
        left_idxs = np.where(directions==1, start_idxs, end_idxs-1)
        right_idxs = np.where(directions==1, end_idxs-1, start_idxs)
        all_directions[offsets[1:-1]] = left_idxs[1:]-right_idxs[:-1]
        all_directions[0] = left_idxs[0]
        return np.cumsum(all_directions)

    def extract_regions(self, regions):
        starts = regions.starts
        ends = regions.ends
        directions = regions.directions
        start_idxs = np.searchsorted(self._indices, starts, side="right")-1
        end_idxs = np.searchsorted(self._indices, ends, side="left")
        offsets = np.insert(np.cumsum(end_idxs-start_idxs), 0, 0)
        slice_indexes=self._get_slice_indexes(start_idxs, end_idxs, directions, offsets)
        values = self._values[slice_indexes]
        all_directions = broadcast(directions, offsets)
        slice_indexes[all_directions==-1] += 1
        indices = np.insert(self._indices, self._indices.size, self._size)[slice_indexes]

        all_starts = broadcast(starts, offsets)
        all_ends = broadcast(ends, offsets)
        transformed_indices = np.where(all_directions==1, indices-all_starts, all_ends-indices)
        transformed_indices[offsets[:-1]]=0
        # assert np.all(np.diff(transformed_indices)), (transformed_indices, regions)
        return BedGraphArray(transformed_indices, values, ends-starts, offsets)

    def threshold(self, value):
        over = self._values >= value
        changes = np.flatnonzero(over[1:] != over[:-1])
        starts = self._indices[changes[::2]+1]
        ends = self._indices[changes[1::2]+1]
        return Regions(starts, ends)

    def _getslice(self, slice_obj):
        assert slice_obj.step is None or slice_obj.step in (1, -1), slice_obj
        start = slice_obj.start or 0
        start_idx = np.searchsorted(self._indices, start, side="right")
        stop = slice_obj.stop
        if stop is None:
            assert self._size is not None
            stop = self._size
        assert self._size is None or stop<=self._size, (slice_obj, self)
        end_idx = np.searchsorted(self._indices, stop, side="left")
        start_value = self._values[start_idx-1]
        indices = np.insert(self._indices[start_idx:end_idx]-start, 0, 0)
        values = np.insert(self._values[start_idx:end_idx], 0, start_value)
        new_obj = self.__class__(indices, values, stop-start)
        if slice_obj.step == -1:
            return new_obj.reverse()
        return new_obj

    def get_start_idx(position):
        return np.searchsorted(self._indices, position, side="right")

    def get_end_idx(position):
        return np.searchsorted(self._indices, position, side="left")

    def scale_x(self, size):
        assert np.issubdtype(self._indices.dtype, np.integer), self._indices
        new_indices = (self._indices*size//self._size)
        assert np.issubdtype(new_indices.dtype, np.integer), (self._indices, size, self._size)
        ds = np.concatenate((np.diff(new_indices)>0, [True]))
        return BedGraph(new_indices[ds], self._values[ds], size)

    @classmethod
    def concatenate(cls, bedgraphs):
        offsets = np.cumsum([0] + [bg._size for bg in bedgraphs[:-1]])
        indices = np.concatenate([bg._indices+offset for (bg, offset) in zip(bedgraphs, offsets)])
        assert np.all(np.diff(indices)>0), (bedgraphs, indices, offsets)
        values = np.concatenate([bg._values for bg in bedgraphs])
        size = sum(bg._size for bg in bedgraphs)
        assert np.all(indices<size), (size, offsets, indices)
        return cls(indices, values, size)

    def update_dense_diffs(self, diffs):
        assert np.issubdtype(self._indices.dtype, np.integer), self._indices
        diffs[0]+=self._values[0]
        diffs[self._indices[1:]]+= np.diff(self._values)

    def __getitem__(self, index):
        if isinstance(index, slice):
            return self._getslice(index)
        if isinstance(index, list):
            return [self._getitem(i) for i in index]
        return self._getitem(index)

    def __eq__(self, other):
        t = np.all(self._indices==other._indices)
        return t and np.all(self._values==other._values)
    
    def __repr__(self):
        return "BG(%s, %s, %s)" % (self._indices, self._values, self._size)

class BedGraphArray:
    def __init__(self, indices, values, sizes, offsets):
        self._indices = np.asanyarray(indices)
        assert np.all(self._indices>=0), self._indices.dtype
        self._values = np.asanyarray(values)
        self._sizes = np.asanyarray(sizes)
        assert np.all(self._sizes>0), sizes

        assert offsets[-1]==self._indices.size, (offsets[-1], self._indices.size)
        assert np.all(np.diff(offsets)>0), offsets
        self._offsets = np.asanyarray(offsets)

    def __repr__(self):
        return f"BGA({self._indices},{self._values}, {self._sizes}, {self._offsets})"

    def __eq__(self, other):
        t = np.all(self._indices==other._indices)
        t &= np.all(self._values==other._values)
        t &= np.all(self._sizes==other._sizes)
        t &= np.all(self._offsets==other._offsets)
        return t

    def scale_x(self, size):
        assert size > 0 
        all_sizes=broadcast(self._sizes, self._offsets)
        new_indices = (self._indices*size//all_sizes)
        mask = np.concatenate((np.diff(new_indices)>0, [True]))
        mask[self._offsets[1:]-1] = True
        counts = np.cumsum(mask)
        new_offsets = np.insert(counts[self._offsets[1:]-1], 0, 0)
        return BedGraphArray(new_indices[mask], self._values[mask], size*np.ones_like(self._sizes), new_offsets)

    def update_dense_diffs(self, diffs, rows):
        assert rows.size == self._offsets.size-1, (rows.size, self._offsets.size-1)
        ncols = diffs.shape[1]
        all_rows = broadcast(rows, self._offsets)
        composite_indexes = all_rows*ncols + self._indices
        args = np.argsort(composite_indexes, kind="mergesort")
        indices = composite_indexes[args]
        index_changes = np.insert(indices[:-1] != indices[1:], indices.size-1, True)
        value_diffs = np.insert(np.diff(self._values), 0, self._values[0])
        value_diffs[self._offsets[:-1]] = self._values[self._offsets[:-1]]
        value_diffs = value_diffs[args]
        totals = np.insert(np.cumsum(value_diffs)[index_changes], 0, 0)
        total_diffs = np.diff(totals)
        used_indices = indices[index_changes]
        diffs[used_indices//ncols, used_indices % ncols] += total_diffs

    def _col_sum(self):
        assert np.all(self._sizes==self._sizes[0]), self._sizes
        args = np.argsort(self._indices, kind="mergesort")
        indices = self._indices[args]
        index_changes = np.insert(indices[:-1] != indices[1:], indices.size-1, True)
        value_diffs = np.insert(np.diff(self._values), 0, self._values[0])[args]
        value_diffs[:self._offsets.size-1] = self._values[self._offsets[:-1]]
        values = np.cumsum(value_diffs)[index_changes]
        indices = indices[index_changes]
        return BedGraph(indices, values, size=self._sizes[0])

    def join_rows(self, offsets):
        cum_sizes = np.insert(np.cumsum(self._sizes), 0, 0)
        new_starts = cum_sizes[:-1]-broadcast(cum_sizes[offsets[:-1]], offsets)
        new_indices = self._indices+broadcast(new_starts, self._offsets)
        new_sizes = np.diff(cum_sizes[offsets])
        new_offsets = self._offsets[offsets]
        return self.__class__(new_indices, self._values, new_sizes, new_offsets)

    def sum(self, axis=None):
        assert axis in (1, None)
        if axis == 1:
            return self._col_sum()

    def extract_regions(self, regions):
        assert regions.starts.size == self._sizes.size
        assert np.all(regions.directions==1), "only positive regions supported"
        all_starts = broadcast(regions.starts, self._offsets)
        all_ends = broadcast(regions.ends, self._offsets)
        mask = all_starts < self._indices
        mask &= self._indices < all_ends
        last_entries = self._offsets[1:]-1
        mask[:-1] |= mask[1:] # Include previous index
        mask[last_entries] |= self._indices[last_entries] < all_starts[last_entries]
        offsets = np.insert(np.cumsum(mask), 0, 0)[self._offsets]
        indices = self._indices[mask]-all_starts[mask]
        assert np.all(indices.size > offsets[:-1]), (offsets, mask, regions, self[0])
        indices[offsets[:-1]] = 0
        values = self._values[mask]
        print(mask, indices, values)
        return self.__class__(indices, values, regions.sizes(), offsets)

    def __getitem__(self, idx):
        assert idx<self._offsets.size-1
        s = slice(self._offsets[idx], self._offsets[idx+1])
        return BedGraph(self._indices[s], self._values[s], self._sizes[idx])

    def __iter__(self):
        return (self[i] for i in range(self._sizes.size))

    @classmethod
    def vstack(cls, arrays):
        index_counts = [a._offsets[-1] for a in arrays]
        offset_offsets = np.cumsum([0] + index_counts)[:-1]
        offset_list = [a._offsets[:-1]+o for a, o in zip(arrays, offset_offsets)]
        last_offset = arrays[-1]._offsets[-1]+offset_offsets[-1]
        new_offsets = np.concatenate(offset_list + [[last_offset]])
        print(new_offsets)
        return cls(np.concatenate([a._indices for a in arrays]),
                   np.concatenate([a._values for a in arrays]),
                   np.concatenate([a._sizes for a in arrays]),
                   new_offsets)

    @classmethod
    def from_bedgraphs(cls, bedgraphs):
        sizes = [bg._size for bg in bedgraphs]
        offsets = np.cumsum([0] + [bg._indices.size for bg in bedgraphs])
        values = np.concatenate([bg._values for bg in bedgraphs])
        indices = np.concatenate([bg._indices for bg in bedgraphs])
        return cls(indices, values, sizes, offsets)
        
    def piecewise_scale(self, break_points, new_sizes):
        new_indices = np.empty_like(self._indices)
        broad_casted = (broadcast(bp, self._offsets) for bp in break_points)
        new_offsets = np.cumsum([0] + list(new_sizes))
        for (all_starts, all_ends), ns, offset in zip(pairwise(broad_casted), new_sizes, new_offsets):
            mask = self._indices>=all_starts
            mask &= self._indices<all_ends
            new_indices[mask] = (offset+(self._indices-all_starts)*ns//(all_ends-all_starts))[mask]

        mask = np.concatenate((np.diff(new_indices)>0, [True]))
        mask[self._offsets[1:]-1] = True
        counts = np.cumsum(mask)
        new_offsets = np.insert(counts[self._offsets[1:]-1], 0, 0)
        return BedGraphArray(new_indices[mask], self._values[mask], sum(new_sizes)*np.ones_like(self._sizes), new_offsets)
