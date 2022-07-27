import Blueprint
const bp = Blueprint

const length = 2.30
const inner_overall_width = 1.25
const inner_board_count = 3

const inner_board_raw_width = inner_overall_width/inner_board_count

const left_wall = bp.plane_at_pos([1.0, 0.0, 0.0], [0.0, 0.0, 0.0])
const right_wall = bp.plane_at_pos([-1.0, 0.0, 0.0], [length, 0.0, 0.0])

const strong_beam = bp.beam_specs(3.0, 4.0)

function inner_insulation_board(offset)
    
end

function render()
    x = inner_insulation_board(0)
    #inner_boards = bp.group(inner_insulation_board(0))
    #return bp.group(inner_boards)
end

render()
