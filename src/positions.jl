struct Position
    seq::String
    pos::Int
end
Position(x::CSV.Row) = Position(string(x.seq), Int(x.pos_start))
