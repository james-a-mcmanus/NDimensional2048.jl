import Base: size, getindex, setindex, show
using SparseArrays
using StatsBase

abstract type AbstractBoard{T,N} <: AbstractArray{T,N} end

mutable struct Board{T} <: AbstractBoard{T,2}
	plates::Matrix{T}
	frees::Array{CartesianIndex{2}}
end

struct SparseBoard{T} <: AbstractBoard{T,2}
	plates::SparseMatrixCSC{T}
end

const StartingFillProportion = 0.5

function SparseBoard(N::Int)
	plates = spzeros(Int, N,N)
	return SparseBoard(plates)
end

function Board(N::Int)
	plate = zeros(Int, N, N)
	free = reduce(vcat, collect(CartesianIndices(plate)))
	b = Board(plate, free)
	start_inds = StatsBase.sample(CartesianIndices(b), floor(Int, N^2 * StartingFillProportion))
	add_number!(b, start_inds)
	return b
end

Base.size(b::AbstractBoard) = size(b.plates)

Base.getindex(b::AbstractBoard, i...) = getindex(b.plates,i...)
Base.setindex!(b::AbstractBoard, i...) = setindex!(b.plates,i...)

frees(b) = iszero.(b)
frees(b::Board) = b.frees
frees(b::SparseBoard) = findall(iszero, b) # slow

plates(b::AbstractBoard) = b.plates

function rand_free_plate(b::Board)
	free_plates = frees(b)
	isempty(free_plates) && return nothing
	return rand(1:length(free_plates))
end

function rand_free_plate(b::SparseBoard; limit=1000)
	counter = 0
	while true
		i = rand(1:length(b))
		counter += 1
		counter > limit && break
		iszero(b[i]) && return i 
	end
	return rand(findall(!iszero, b))
end

function next_turn!(board)
	plate_ind = rand_free_plate(board)
	isnothing(plate_ind) && return 1
	add_number!(board, plate_ind)
	return 0
end

function add_number!(board::Board, ind::Int)
	plates(board)[frees(board)[ind]] = new_number(board)
	deleteat!(frees(board), ind)
end

function add_number!(board::Board, ind::CartesianIndex)
	plates(board)[ind] = new_number(board)
	freeind = findfirst(x->x==ind, frees(board))
	!isnothing(freeind) && deleteat!(frees(board), freeind)
end

function add_number!(board::SparseBoard, ind)
	plates(board)[ind] = new_number(board)
end
function add_number!(b::Board, ind::AbstractArray)
	plates(b)[ind] .= new_number(b)
end

new_number(a) = 2

function move!(board, iterover, rev)
	board_length = size(board,1)
	for slice in iterover(board)
		squeeze_slice!(slice, rev)
		calc_collisions!(slice, rev ? same : reverse)
		squeeze_slice!(slice, rev)
	end	
end

function move!(board::SparseBoard, iterover, rev)

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

function i_to_cart(a, i)
	inds = CartesianIndices(a)
	return inds[i]
end