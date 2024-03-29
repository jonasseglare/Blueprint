using Revise
import Blueprint
using Setfield

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

const fitness_margin = 0.001
const left_wall_marg = bp.translate_normalized(left_wall, fitness_margin)
const right_wall_marg = bp.translate_normalized(right_wall, fitness_margin)

const base = bp.plane_at_pos([0.0, 0.0, 1.0], [0.0, 0.0, 0.0])

const crop_margin = 0.1

const left_cut = bp.NamedPlane(:left_cut, bp.translate_normalized(left_wall, -crop_margin))
const right_cut = bp.NamedPlane(:right_cut, bp.translate_normalized(right_wall, -crop_margin))

const strong_specs = bp.set_label(bp.beam_specs(0.03, 0.04), "Strong beam")
const square_specs = bp.set_label(bp.beam_specs(0.03, 0.03), "Square beam")
const plywood_sheet_specs = bp.set_label(bp.set_cutting_plan_key(bp.beam_specs(1.0, 0.003), :beam_Y_lower), "Plywood sheet")

# Cut the ends of an object resting on the shelf
function shelf_cut(x)
    return bp.cut(left_cut, bp.cut(right_cut, x))
end

function supporting_beam(dir, specs=strong_specs)
    return bp.push_against(base, bp.orient_beam(bp.new_beam(specs), dir, bp.local_y_dir))
end

function inner_insulation_board(offset)
    test_plane = bp.NamedPlane(:close_plane, bp.plane_at_pos([0.0, 1.0, 0.0], [0.0, 0.0, 0.0]))
    close_plane = bp.NamedPlane(:close_plane, bp.plane_at_pos([0.0, 1.0, 0.0], [0.0, offset + fitness_margin, 0.0]))
    far_plane = bp.NamedPlane(:far_plane, bp.plane_at_pos([0.0, -1.0, 0.0], [0.0, offset + inner_board_raw_width - fitness_margin, 0.0]))

    function cut_close_far(x)
        return bp.cut_many([close_plane, far_plane], x)
    end

    close_beam = shelf_cut(bp.push_against(close_plane, supporting_beam([1.0, 0.0, 0.0])))
    far_beam = shelf_cut(bp.push_against(far_plane, supporting_beam([1.0, 0.0, 0.0])))

    close_inner = bp.get_tangent_cutting_plane(:close_inner, close_beam, [0.0, 1.0, 0.0])
    far_inner = bp.get_tangent_cutting_plane(:far_inner, far_beam, [0.0, -1.0, 0.0])

    support_specs = bp.set_diagram_strategies(strong_specs, bp.default_diagram_strategy, bp.basic_diagram_strategy(:individual_first, :close_inner, bp.x_vector_3d))
    
    connecting_beams = bp.with_memberships(bp.cut_many([close_inner, far_inner], bp.beam_array_between_planes(left_cut.plane, right_cut.plane, supporting_beam([0.0, 1.0, 0.0], support_specs), 4)), :connecting_beams)

    x_drill_planes = bp.spaced_planes(bp.pinch_beams(connecting_beams, bp.x_vector_3d), 1)
    z_drill_planes = bp.spaced_planes(bp.pinch_component(connecting_beams, bp.z_vector_3d), 2, 0.7)

    drill_radius = 0.003
    common_drill_specs = bp.DrillSpecs(drill_radius, bp.LabelSpec("Common drill", "C"), true, true)
    
    close_conbeam_drills = bp.generate_drills(bp.y_vector_3d, x_drill_planes, z_drill_planes, common_drill_specs)
    
    connecting_beams = bp.drill(connecting_beams, close_conbeam_drills)

    plywood_sheet_proto = shelf_cut(cut_close_far(bp.orient_beam(bp.new_beam(plywood_sheet_specs), [1.0, 0.0, 0.0], bp.local_y_dir)))
    
    down = [0.0, 0.0, -1.0]
    up = -down
    plywood_sheet_above = bp.push_component_against_component(close_beam, down, plywood_sheet_proto)
    plywood_sheet_under = bp.push_component_against_component(close_beam, up, plywood_sheet_proto)

    alignment_beam_proto = cut_close_far(bp.push_component_against_component(plywood_sheet_under, up, bp.orient_beam(bp.new_beam(square_specs), [0.0, 1.0, 0.0], bp.local_y_dir)))
    left_alignment_beam = bp.push_against(left_wall_marg, alignment_beam_proto)
    right_alignment_beam = bp.push_against(right_wall_marg, alignment_beam_proto)

    alignment_beams = bp.membership_group(:alignment, [left_alignment_beam, right_alignment_beam])
    
    covers = bp.membership_group(:covers, [plywood_sheet_above, plywood_sheet_under])
    support_beams = bp.drill(bp.membership_group(:main_support, [close_beam, far_beam]), close_conbeam_drills)
    
    structure = bp.membership_group(:structure, [support_beams, connecting_beams, alignment_beams])
    result = bp.group([covers, structure])
    return result
end

function render()
    full_design = bp.group(inner_insulation_board(0))

    vertical_view = bp.ProjectedView("View from above", bp.flip(base), [1.0, 0.0, 0.0], nothing, bp.member_of_union(:structure))
    side_view = bp.ProjectedView("Side view", bp.plane_at_pos([0.0, 1.0, 0.0], [0.0, 0.0, 0.0]), [1.0, 0.0, 0.0], nothing, bp.member_of_union(:structure))
    
    report = bp.basic_report("Insulation boards", full_design)
    bp.push_sub_model!(report, bp.SubModel("Structure", "structure.stl", (memberships) -> :structure in memberships))
    bp.push_view!(report, vertical_view)
    bp.push_view!(report, side_view)
    
    doc = bp.make("output/insulation_boards", report)
    bp.render_markdown(doc)
    bp.render_html(doc)
    print("Rendered it")
end

render()
