using Revise
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
const plywood_sheet_specs = bp.beam_specs(1.0, 0.003)

function shelf_cut(x)
    return bp.cut(left_cut, bp.cut(right_cut, x))
end

function supporting_beam(dir)
    return bp.push_against(base, bp.orient_beam(bp.new_beam(strong_specs), dir, bp.local_y_dir))
end




function beam_array_between_planes(start_plane, end_plane, beam_prototype, count)
    @assert 2 <= count
    start_plane = bp.normalize_plane(start_plane)
    end_plane = bp.normalize_plane(end_plane)
    start_translation = bp.push_against_transform(start_plane, beam_prototype).translation
    end_translation = bp.push_against_transform(end_plane, beam_prototype).translation
    step_count = count - 1
    step = (1.0/step_count)*(end_translation - start_translation)
    return bp.group([bp.transform(bp.rigid_transform_from_translation(start_translation + step*i), beam_prototype) for i in 0:step_count])
end

function inner_insulation_board(offset)
    close_plane = bp.NamedPlane(:close_plane, bp.plane_at_pos([0.0, 1.0, 0.0], [0.0, offset, 0.0]))
    far_plane = bp.NamedPlane(:far_plane, bp.plane_at_pos([0.0, -1.0, 0.0], [0.0, offset + inner_board_raw_width, 0.0]))

    function cut_close_far(x)
        return bp.cut(close_plane, bp.cut(far_plane, x))
    end

    close_beam = shelf_cut(bp.push_against(close_plane, supporting_beam([1.0, 0.0, 0.0])))
    far_beam = shelf_cut(bp.push_against(far_plane, supporting_beam([1.0, 0.0, 0.0])))

    connecting_beams = cut_close_far(beam_array_between_planes(left_cut.plane, right_cut.plane, supporting_beam([0.0, 1.0, 0.0]), 4))

    plywood_sheet = shelf_cut(cut_close_far(bp.orient_beam(bp.new_beam(plywood_sheet_specs), [1.0, 0.0, 0.0], bp.local_y_dir)))
    
    return bp.group([#close_beam, far_beam, connecting_beams,
                     bp.push_component_against_component(close_beam, [0.0, 0.0, -1.0], plywood_sheet)])
end

function render()
    full_design = bp.group(inner_insulation_board(0))

    report = bp.basic_report("Insulation boards", full_design)
    doc = bp.make("output/insulation_boards", report)
    bp.render_markdown(doc)
    bp.render_html(doc)
end

render()
