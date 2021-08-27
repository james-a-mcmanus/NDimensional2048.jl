import Base: size, getindex, setindex

mutable struct Board{T} <: AbstractArray{T,2}
	plates::Matrix{T}
	frees::Array{CartesianIndex{2}}
end

function Board(N::Int)
	plates = zeros(Int, N, N)
	frees = reduce(vcat, collect(CartesianIndices(plates)))
	return Board(plates, frees)
end

Base.size(b::Board) = size(b.plates)

Base.getindex(b::Board, i...) = getindex(b.plates,i...)
Base.setindex!(b::Board, i...) = setindex!(b.plates,i...)

frees(b) = iszero.(b)
frees(b::Board) = b.frees
plates(b::Board) = b.plates


function next_turn!(board)
	free_plates = frees(board)
	isempty(free_plates) && return 1
	i = rand(1:length(free_plates))
	chosen_plate = free_plates[i]
	add_number!(board, chosen_plate)
	deleteat!(free_plates, i)
	return 0
end

function add_number!(board, ind)
	plates(board)[ind] = 2
end

function move!(board, iterover, rev)
	board_length = size(board,1)
	for slice in iterover(board)
		squeeze_slice!(slice, rev)
		calc_collisions!(slice, rev ? same : reverse)
		squeeze_slice!(slice, rev)
	end	
end

up!(board) = move!(board, eachcol, true)
down!(board) = move!(board, eachcol, false)
right!(board) = move!(board, eachrow, false)
left!(board) = move!(board, eachrow, true)

front(itr, n=1) = Iterators.take(itr, length(itr)-n)
back(itr, n=1) = Iterators.drop(itr, n)

same(a) = a
same(a,b) = (a,b)

"""
find any that are next to each other and combine them.
"""
function calc_collisions!(sl, rev)

	N = length(sl)
	r0 = rev(2:N)
	r1 = rev(1:N-1)
	(r1, r0) = rev((r0, r1))

	for i in eachindex(r0, r1)
		sl[r0[i]] == 0 && continue
		sl[r0[i]] != sl[r1[i]] && continue
		sl[r0[i]] += sl[r1[i]]
		sl[r1[i]] = 0
	end
end

"""
move all the filled pieces along.
"""
function squeeze_slice!(slice, rev)
	slice .= slice[sortperm(slice .!= 0, rev=rev)]
end

"""
get the score for the board
"""
score(board) = sum(board)


"""
get a character from the keyboard
"""
function getc1()
   ret = ccall(:jl_tty_set_mode, Int32, (Ptr{Cvoid},Int32), stdin.handle, true)
   ret == 0 || error("unable to switch to raw mode")
   c = read(stdin, Char)
   ccall(:jl_tty_set_mode, Int32, (Ptr{Cvoid},Int32), stdin.handle, false)
   c
end


"""
play a game.
"""
function play(;N=1)
	b = Board(N)
	show_board(b)
	while true
		opt = getc1()
		is_arrow(opt) || continue
		opt == 'q' && break
		make_move!(b, opt) || break
		next_turn!(b) || break
		show_board(b)
	end
end


function make_move(board, key)
	key == '\e'
		left!(board)
	#key = 
end