from __future__ import annotations

from typing import Iterator, List, Type



class ChromoRegionMergeException(Exception):
	pass

class ChromoRegion:
	__slots__ = ["chromo", "start", "_end", "_length"]

	"""Describes a Chromosome region. Index is 0-based and half-closed"""
	def __init__(self, chromo: str, start: int, end: int) -> None:
		assert start <= end
		self.chromo = chromo
		self.start = start
		self._end = end
		self._length = end - start

	@property
	def end(self) -> int:
		return self._end

	@end.setter
	def end(self, value: int) -> None:
		assert value >= self.start

		self._end = value
		self._length = self._end - self.start

	def contiguousWith(self, o: ChromoRegion) -> bool:
		if self.chromo != o.chromo:
			return False

		if self.start < o.start:
			return self.end >= o.start
		else:
			return o.end >= self.start

	def __add__(self, o: ChromoRegion) -> ChromoRegion:
		if not self.contiguousWith(o):
			raise ChromoRegionMergeException(f"Regions are not contiguous: {self.chromo}:{self.start}-{self.end}, {o.chromo}:{o.start}-{o.end}")

		return ChromoRegion(self.chromo, min(self.start, o.start), max(self.end, o.end))

	def __sub__(self, o: ChromoRegion) -> List[ChromoRegion]:
		if not self.contiguousWith(o):
			return [ChromoRegion(self.chromo, self.start, self.end)]

		if o.start <= self.start and o.end >= self.end:
			return []
		elif o.start > self.start and o.end < self.end:
			return [ChromoRegion(self.chromo, self.start, o.start), ChromoRegion(self.chromo, o.end, self.end)]
		elif o.start <= self.start and o.end < self.end:
			return [ChromoRegion(self.chromo, o.end, self.end)]
		elif o.start > self.start and o.end >= self.end:
			return [ChromoRegion(self.chromo, self.start, o.start)]

		# Should not happen, the above conditions are exhaustive
		return []

	def __len__(self) -> int:
		return self._length

	def __eq__(self, o: object) -> bool:
		if not isinstance(o, ChromoRegion):
			return NotImplemented

		return (self.chromo == o.chromo
			and self.start == o.start
			and self.end == o.end
			and self._length == o._length)

	# __lt__ has been implemented specifically for sorting -- both the "list.sort()" method
	# and "sorted" function use __lt__ for comparing objects.
	def __lt__(self, o: object) -> bool:
		if not isinstance(o, ChromoRegion):
			return NotImplemented

		if self.chromo == o.chromo:
			return self.start < o.start

		selfChromo = self.chromo
		if selfChromo.startswith("chr"):
			selfChromo = selfChromo[3:]

		otherChromo = o.chromo
		if otherChromo.startswith("chr"):
			otherChromo = otherChromo[3:]

		if selfChromo.isdigit() and otherChromo.isdigit():
			return int(selfChromo) < int(otherChromo)

		if not selfChromo.isdigit() and otherChromo.isdigit():
			return False # sort letters higher than numbers

		if selfChromo.isdigit() and not otherChromo.isdigit():
			return True # sort letters higher than numbers

		return selfChromo < otherChromo # sort lexigraphically

	def __repr__(self) -> str:
		return f"({self.chromo}:{self.start}-{self.end})"


# Not quite a set in the mathematical sense.
class ChromoRegionSet:
	__slots__ = ["regions", "_chromoSet", "_chromos", "cumulativeRegionSize", "_chromoOrderDirty"]

	regions: list[ChromoRegion]
	_chromoSet: set[str]
	_chromos: list[str]
	cumulativeRegionSize: int
	_chromoOrderDirty: bool

	def __init__(self, regions: List[ChromoRegion]=None) -> None:
		self.regions = []
		self._chromoSet = set()
		self._chromos = []
		self.cumulativeRegionSize = 0
		self._chromoOrderDirty = False
		if regions is not None:
			self.regions = regions
			for region in self.regions:
				if region.chromo not in self._chromoSet:
					self._chromoSet.add(region.chromo)
					self._chromos.append(region.chromo)
				self.cumulativeRegionSize += len(region)

	def addRegion(self, region: ChromoRegion) -> None:
		self._chromoSet.add(region.chromo)
		self.regions.append(region)
		self.cumulativeRegionSize += len(region)
		self._chromoOrderDirty = True

	@property
	def chromos(self):
		if self._chromoOrderDirty:
			chromoSet = set()
			chromos = []
			for region in self.regions:
				if region.chromo not in chromoSet:
					chromoSet.add(region.chromo)
					chromos.append(region.chromo)
			self._chromos = chromos
		self._chromoOrderDirty = False
		return self._chromos

	def sortRegions(self) -> None:
		self.regions.sort()
		self._chromoOrderDirty = True

	def mergeRegions(self) -> None:
		if len(self.regions) == 1:
			return

		self.sortRegions()

		mergedRegions = []
		currentRegion = self.regions[0]

		for region in self.regions[1:]:
			if currentRegion.contiguousWith(region):
				currentRegion = currentRegion + region
			else:
				mergedRegions.append(currentRegion)
				currentRegion = region
		mergedRegions.append(currentRegion)

		mergedCumulativeLength = 0
		for region in mergedRegions:
			mergedCumulativeLength += len(region)

		self.regions = mergedRegions
		self.cumulativeRegionSize = mergedCumulativeLength

	def __add__(self, o: ChromoRegionSet) -> ChromoRegionSet:
		"""Simple conacatenation of sets. No region merging is done."""
		newRegionSet = ChromoRegionSet(self.regions + o.regions)
		newRegionSet.sortRegions()

		return newRegionSet

	def __sub__(self, o: ChromoRegionSet) -> ChromoRegionSet:
		regionWorkingSet = self.regions
		for removeRegion in o:
			tempRegions = []
			for region in regionWorkingSet:
				tempRegions.extend(region - removeRegion)
			regionWorkingSet = tempRegions

		return ChromoRegionSet(regionWorkingSet)

	def __len__(self) -> int:
		return len(self.regions)

	def __iter__(self) ->  Iterator[ChromoRegion]:
		for region in self.regions:
			yield region

	def __eq__(self, o: object) -> bool:
		if not isinstance(o, ChromoRegionSet):
			return NotImplemented

		if len(self.regions) != len(o.regions):
			return False

		if self.cumulativeRegionSize != o.cumulativeRegionSize:
			return False

		for selfRegion, oRegion in zip(sorted(self), sorted(o)):
			if selfRegion != oRegion:
				return False

		return True

	def __repr__(self) -> str:
		return f"[{', '.join([str(region) for region in self.regions])}]"

	@classmethod
	def loadBed(cls: Type[ChromoRegionSet], filename: str) -> ChromoRegionSet:
		regionSet = ChromoRegionSet()

		with open(filename) as regionFile:
			regionLines = regionFile.readlines()

		for line in regionLines:
			temp = line.split()
			regionSet.addRegion(ChromoRegion(temp[0], int(temp[1]), int(temp[2])))

		return regionSet
