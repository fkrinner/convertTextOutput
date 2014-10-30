class vector:
	def __init__(self, vec):
		self.vec = vec

	def __add__(self,other):
		if len(self) != len(other):
			print "Can't add vectors of different length"
			raise IndexError
		return vector([self[i]+other[i] for i in range(len(self))])

	def __sub__(self,other):
		if len(self) != len(other):
			print "Can't subtract vectors of different length"
			raise IndexError
		return vector([self[i]-other[i] for i in range(len(self))])

	
	def __mul__(self,other):
		if isinstance(other, vector):
			if len(self) != len(other):
				print "Can't multiply vectors of different length"
				raise IndexError
			sp = 0.
			for i in range(len(self)):
				sp+=self[i]*other[i]
			return sp
		else:
			return vector([other * self[i] for i in range(len(self))])

	def __div__(self,other):
		return self.__mul__(1./other)

	def __rmul__(self,other):
		return self.__mul__(other)


	def __abs__(self):
		return (self.conjugate()*self)**.5

	def __str__(self):
		return "vector"+str(self.vec)
	
	def __getitem__(self,i):
		return self.vec[i]
		

	def __len__(self):
		return len(self.vec)

	def __iter__(self):
		for val in self.vec:
			yield val

	def normalize(self):
		self/=abs(self)

	def conjugate(self):
		return vector([val.conjugate() for val in self.vec])


def remove_vector_from_onb(onb,rem):
	size = len(rem)
	for vec in onb:
		if not len(vec) == size:
			print "length of vectors not equal"
			raise IndexError
	norm_rem = rem/abs(rem)
	onb_norm = []
	for vec in onb:
		onb_norm.append(vec/abs(vec))
	ret_basis = []
	for vec in onb_norm:
		new_basis_vector = vec - norm_rem*(norm_rem*vec)
#		print "after rem norm_rem:",norm_rem*new_basis_vector,new_basis_vector
		for prev in ret_basis:
			new_basis_vector = new_basis_vector - prev*(prev*new_basis_vector)
#			print "after rem_prev:",norm_rem*new_basis_vector,new_basis_vector
		try:
			newabs = abs(new_basis_vector)
			if newabs < 1.E-14:
				raise ZeroDivisionError
			new_basis_vector/=newabs
#			print "append:",new_basis_vector
			ret_basis.append(new_basis_vector)
		except ZeroDivisionError:
			pass
	return ret_basis




if __name__ == "__main__":
	x = vector([1.,0.,0.])
	y = vector([0.,1.,0.])
	z = vector([0.,0.,1.])

	bas = [x,y,z]

	remvec = vector([2.**-.5,2**-.5,0.])

	nb = remove_vector_from_onb(bas,remvec)
	for vec in nb:
		print vec,abs(vec)
	print nb[0]*nb[1]






