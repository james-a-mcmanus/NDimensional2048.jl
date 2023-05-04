import Base: size, getindex, setindex, show
using SparseArrays
using StatsBase
using Crayons
using Colors

abstract type AbstractBoard{T,N} <: AbstractArray{T,N} end

mutable struct Board{T} <: AbstractBoard{T,2}
	plates::Matrix{T}
	frees::Array{CartesianIndex{2}}
end

struct BoardND{T,N} <: AbstractBoard{T,N}
	plates::Array{T,N}
end

struct SparseBoard{T} <: AbstractBoard{T,2}
	plates::SparseMatrixCSC{T}
end

const StartingFillProportion = 0.5
const panel_colors = Dict{Number, Crayon}()
function default_colors!(cd)

	maxnum = 12
	startcolor = HSL(10,.8,.5)
	endcolor=HSL(360,.8,.5)
	colorange = RGB.(range(startcolor, endcolor, length=maxnum))
	for i in 1:maxnum
		c = colorange[i]
		cd[2^i] = Crayon(background=(round(Int,red(c)*255), round(Int,green(c)*255), round(Int,blue(c)*255)), bold=true, reset=true)
	end
	cd[NaN] = Crayon(foreground=(100,100,100), background=(100,100,100), bold=true, reset=true)
	cd[0] = Crayon(foreground=(100,100,100), background=(100,100,100), bold=true, reset=true)
end
default_colors!(panel_colors)


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

function BoardND(W::Int, N::Int)
	sizeof(Int) * W^N > Sys.free_memory() && error("Requested size too large, try SparseBoardND instead.")
	plate = zeros(Int, repeat([W],N)...)
	return BoardND(plate)
end

Base.size(b::AbstractBoard) = size(b.plates)

Base.getindex(b::AbstractBoard, i...) = getindex(b.plates,i...)
Base.setindex!(b::AbstractBoard, i...) = setindex!(b.plates,i...)

show(x::SparseBoard) = show(plates(x))
show(io::IO, ::MIME"text/plain", x::SparseBoard) = show(io, plates(x))

frees(b) = iszero.(b)
frees(b::Board) = b.frees
frees(b::SparseBoard) = findall(iszero, b) # slow
frees(b::BoardND) = findall(iszero, b)

plates(b::AbstractBoard) = b.plates

function rand_free_plate(b::Board)
	free_plates = frees(b)
	isempty(free_plates) && return nothing
	return rand(1:length(free_plates))
end

function rand_free_plate(b::Union{SparseBoard, BoardND}; limit=1000)
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
function next_turn!(b,N)
	for _ in 1:N
		next_turn!(b)
	end
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

function add_number!(board::Union{SparseBoard, BoardND}, ind)
	plates(board)[ind] = new_number(board)
end
function add_number!(b::Board, ind::AbstractArray)
	plates(b)[ind] .= new_number(b)
end

new_number(a) = 2

"""Shift the pieces along a dimension, merging similar pieces that move into eachother."""
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
function move!(board, direction::AbstractVector)
	abs(sum(direction)) == 1 || error("Not an orthogonal vector!")
	directionispos = sum(direction) > 0
	iterover = CartesianIndices(Tuple(direction[i] ==0 ? s : 1 for (i,s) in enumerate(size(board))))
	for i in iterover
		slice = view(board, _vectorrange(i, abs.(direction), size(board))...)
		squeeze_slice!(slice, directionispos)
		calc_collisions!(slice, directionispos ? same : reverse)
		squeeze_slice!(slice, directionispos)
	end
end

function _vectorrange(start, direction, arraysize)
	index = Vector{Union{Int, UnitRange}}(undef, length(direction))
	for d in 1:length(direction)
		index[d] = direction[d] !=0 ? UnitRange(start[d], start[d] + (direction[d] * arraysize[d] - 1)) : start[d]
	end
	return index
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

function play!(b)
	game_finished = 0
	while game_finished == 0
		render(b)
		mv = readline()
		make_move!(b, mv)
		game_finished = next_turn!(b)
		print("\u1b[1J")
		render(b)
	end
end


function make_move!(board, key)
	key == "\e[D" && left!(board)
	key == "\e[C" && right!(board)
end

function i_to_cart(a, i)
	inds = CartesianIndices(a)
	return inds[i]
end


"""
Display a board
"""
function render(b::Board)

	print_border(b)
	term_size = displaysize(stdout)
	panel_size = _get_panel_size(term_size, size(b))
	inds = eachindex(b)
	for panel in eachrow(inds)
		print(panel_colors[NaN], '|')
		panel_char = _get_panel_char.(panel_size, b[panel])
		for p in panel_char
			print(p[1], p[2])
		end
		print(panel_colors[NaN] ,'|')
		println('\n' * " "^(size(b,2)*4+2))
	end
	print_border(b)
end

function print_border(b)
	print(panel_colors[NaN], "_" ^ (size(b,2)*4+2))
	print('\n')
end

function _check_eol(ind, sz)
	ind[2] == sz[2]
end

function _get_panel_char(ps, n)
	cl = panel_colors[n]
	n = rpad(string(n),4)
	return cl, n
end

function _get_panel_size(ts, bs)
	return 1
end

function testrecall()
	println(",")
	println(",")
end

function clearconsole()
	run(`clear`);
end