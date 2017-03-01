program mesh
implicit none
character(20) name_file
real(8),allocatable:: node(:,:)
integer(4),allocatable:: element(:,:),mass_node(:,:),contact_elem(:,:)
integer(4),allocatable:: contact(:,:),cont_numer(:,:),bound(:,:)
integer(4) i,j,k,l,nn,s,ind,flag(10),f(2),d,d1,d2,n,ll,ii,jj,kk
integer(4) N_node,N_elem,num
integer(4) n_point,n_cell,n_bound

read(*,*) name_file
open(1,file = name_file) 
read(1,*) 
read(1,*) 
read(1,*) 
read(1,*) 
read(1,*) N_node
allocate(node(n_node,4))
node = 0 
do i = 1,N_node
read(1,*) node(i,:)
enddo

read(1,*) 
read(1,*) 
read(1,*) N_elem
allocate(element(N_elem,9))
element = 0

do i = 1,n_elem
	read(1,*) element(i,1:2)
enddo
close(1)

open(1,file=name_file) 
do i = 1, n_node+8
read(1,*) 
enddo 

n_point = 0 
n_cell = 0
n_bound = 0

do i = 1, n_elem
if (element(i,2) == 15) then
	read(1,*) element(i,1:6)
	n_point = n_point+1	
endif
if (element(i,2) == 1) then
	read(1,*) element(i,1:7)
	n_bound = n_bound+1
endif
if (element(i,2) == 3) then
	read(1,*) element(i,1:9)
	n_cell = n_cell+1
endif	
enddo
close(1)

allocate(mass_node(n_node,9))
mass_node = 0 

do i = 1, n_elem
	if (element(i,2) == 3) then
		do j = 6,9
			k = mass_node(element(i,j),1) 
			k = k + 1
			mass_node(element(i,j),1) = k
			mass_node(element(i,j),k+1) = element(i,1)
		enddo
	endif
enddo

allocate(contact_elem(n_elem,30))

contact_elem = 0
do i = 1,n_node
	num = mass_node(i,1) 
	do j = 1,num
		k = mass_node(i,j+1)
		ext:do l = 1,size(contact_elem(k,:))
			if (contact_elem(k,l) == 0) then
				contact_elem(k,l:l+num+1) = mass_node(i,2:num+1)
				exit ext
			endif
		enddo ext
	enddo
enddo

allocate(contact(n_elem,9),cont_numer(n_elem,8))
contact = 0

do i = 1,n_elem
	s = 1 !Начальная позиция первого контакта
	do j = 1,size(contact_elem(1,:))
	
	    nn = 0	!Количество вхождений
		if (i.ne.contact_elem(i,j)) then !Выктдываем сам элемент, оставляем только контактные элементы
			extt: do l = 1,size(contact_elem(1,:)) !Обход по каждому элементу контакта
				if (contact_elem(i,j).ne.0) then !Выкидываем нули
					
					if (contact_elem(i,j) == contact_elem(i,l)) then !Если cсовпали, то +
						nn = nn+1 !Если у нас совпадение, то nn++
					endif
						
					if (nn == 2) then !Если два совпадения, то добавляем контактный элемент
						contact(i,s) = contact_elem(i,l)
						contact_elem(i,l) = 0
						contact(i,9) = s !Sum contact elements
						s = s+1
						exit extt
					endif
				endif
			enddo extt
		endif
	enddo
enddo

do i = 1,n_elem 
	s = 1
	n = 5
	do j = 1,4
		if (contact(i,j)>0) then
		ll = 1
		do k = 6,9
			do l = 6,9
				if (element(i,k) == element(contact(i,j),l)) then
					cont_numer(i,s) = k-5
					s = s+1
					f(ll) = k-5
					ll = ll+1
				endif
			enddo
		enddo
		call posit(f,d)
		contact(i,n) = d
		n = n+1
		endif
	enddo
enddo

allocate(bound(n_bound,7))
k = 1
do i = 1,n_elem
	if (element(i,2) == 1) then
		bound(k,1:3) = element(i,5:7)
		k = k+1
	endif
enddo



do i = 1,n_bound
	ex:do j = n_point+n_bound+1,n_elem
		ll = 1
		do ii = 2,3
			do jj = 6,9
					if (bound(i,ii) == element(j,jj)) then
						f(ll) = jj - 5
						ll = ll+1
					endif	
					if (ll == 3) then
						bound(i,4) = j
						
						call posit(f,d)
						call posit_cont(d,d1)
						bound(i,2) = d
						bound(i,5) = d1
						call posit_cont(d,d2)
						bound(i,7) = d2
							
						do k = 5,contact(j,9)+4
							if (contact(j,k) == d2) bound(i,6) = contact(j,k-4) 
						enddo
						
						exit ex
					endif
				enddo
			enddo 
	enddo ex
enddo					
!				


open(1,file = 'mesh.inp')
write(1,*) n_node
do i = 1,n_node
	write(1,*) node(i,2:)
enddo

write(1,*) n_bound
do i = 1,n_bound
	write(1,*) bound(i,1:2),bound(i,4)-n_point-n_bound,bound(i,5),bound(i,6)-n_point-n_bound,bound(i,7)
enddo	

write(1,*) n_cell
do i = n_point+n_bound+1,n_cell+n_point+n_bound
	write(1,*) element(i,5:9)
enddo

write(1,*)  
do i = n_point+n_bound+1,n_cell+n_point+n_bound
	do j = 1,4
        if (contact(i,j).ne.0) contact(i,j) = contact(i,j)-n_point-n_bound
	enddo	
	write(1,*) contact(i,:)  
enddo

contains 

subroutine posit(f,d)
integer(4),intent(in):: f(2)
integer(4),intent(out):: d 
if (((f(1) == 1) .and. (f(2) == 2)).or.((f(1) == 2) .and. (f(2) == 1))) d = 1
if (((f(1) == 2) .and. (f(2) == 3)).or.((f(1) == 3) .and. (f(2) == 2))) d = 2
if (((f(1) == 3) .and. (f(2) == 4)).or.((f(1) == 4) .and. (f(2) == 3))) d = 3
if (((f(1) == 1) .and. (f(2) == 4)).or.((f(1) == 4) .and. (f(2) == 1))) d = 4

end subroutine posit

subroutine posit_cont(d1,d2)
integer(4),intent(in):: d1
integer(4),intent(out):: d2
if (d1 == 1) d2 = 3
if (d1 == 2) d2 = 4
if (d1 == 3) d2 = 1
if (d1 == 4) d2 = 2
end subroutine posit_cont

end  
