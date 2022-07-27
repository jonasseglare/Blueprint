import Blueprint
const bp = Blueprint


#################################
# Coordinate system seen from above when entering
#
# Y
# ^
# |
# |
# o----> X
#
#
#  ^
#  | Entering from here.
#

const length = 2.30
const inner_overall_width = 1.25
const inner_board_count = 3

const inner_board_raw_width = inner_overall_width/inner_board_count

# Left wall plane when entering the cellar
const left_wall = bp.plane_at_pos([1.0, 0.0, 0.0], [0.0, 0.0, 0.0])

# Right wall plane when entering the cellar
const right_wall = bp.plane_at_pos([-1.0, 0.0, 0.0], [length, 0.0, 0.0])
const base = bp.plane_at_pos([0.0, 0.0, 1.0], [0.0, 0.0, 0.0])

const crop_margin = 0.1

const left_cut = bp.NamedPlane(:left_cut, bp.translate_normalized(left_wall, -crop_margin))
const right_cut = bp.NamedPlane(:right_cut, bp.translate_normalized(right_wall, -crop_margin))

const strong_specs = bp.beam_specs(0.03, 0.04)

function shelf_cut(x)
    return bp.cut(left_cut, bp.cut(right_cut, x))
end

function supporting_beam(dir)
    return bp.push_against(base, bp.orient_beam(bp.new_beam(strong_specs), dir, bp.local_y_dir))
end

function inner_insulation_board(offset)
    close_plane = bp.plane_at_pos([0.0, 1.0, 0.0], [0.0, offset, 0.0])
    far_plane = bp.plane_at_pos([0.0, -1.0, 0.0], [0.0, offset + inner_board_raw_width, 0.0])

    close_beam = shelf_cut(bp.push_against(close_plane, supporting_beam([1.0, 0.0, 0.0])))
    far_beam = shelf_cut(bp.push_against(far_plane, supporting_beam([1.0, 0.0, 0.0])))
    
    return bp.group([close_beam, far_beam])
end

function render()
    full_design = bp.group(inner_insulation_board(0))

    report = bp.basic_report("Insulation boards", full_design)
    doc = bp.make("output/insulation_boards", report)
    bp.render_markdown(doc)
    bp.render_html(doc)
end

render()
