#!/usr/bin/env python3
"""Render a HOOMD GSD trajectory with OVITO."""

from pathlib import Path
import argparse

from ovito.io import import_file
from ovito.vis import TachyonRenderer, Viewport


def render_frame(
    input_file: Path,
    output_file: Path,
    frame: int,
    size: tuple[int, int],
) -> None:
    pipeline = import_file(str(input_file), multiple_frames=True)
    pipeline.add_to_scene()

    # GSD bond topology is preserved by the importer. This modifier only
    # creates visible bonds when a trajectory does not contain topology.
    if pipeline.source.num_frames == 0:
        raise RuntimeError(f"No frames found in {input_file}")

    viewport = Viewport(
        type=Viewport.Type.Perspective,
        camera_dir=(1.5, 1.5, 1.2),
        camera_pos=(0, 0, 0),
    )
    viewport.zoom_all()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    viewport.render_image(
        filename=str(output_file),
        size=size,
        frame=frame,
        renderer=TachyonRenderer(shadows=True),
    )
    pipeline.remove_from_scene()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path, help="Input trajectory.gsd")
    parser.add_argument("-o", "--output", type=Path, default=Path("frame.png"))
    parser.add_argument("--frame", type=int, default=0)
    parser.add_argument("--width", type=int, default=1400)
    parser.add_argument("--height", type=int, default=1000)
    args = parser.parse_args()
    render_frame(
        args.input,
        args.output,
        args.frame,
        (args.width, args.height),
    )
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
