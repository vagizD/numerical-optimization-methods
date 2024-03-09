#pragma once

enum class DiffMethod { Stencil3, Stencil3Extra, Stencil5, Stencil5Extra, FwdADD };

enum class Derivative { X, Y, XX, YY, XY };

enum class Variable { X, Y };
