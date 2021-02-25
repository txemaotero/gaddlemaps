//    Gaddlemaps python module.
//    Copyright (C) 2019-2021 José Manuel Otero Mato, Hadrián Montes Campos, Luis Miguel Varela Cabo
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Affero General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Affero General Public License for more details.
//
//    You should have received a copy of the GNU Affero General Public License
//    along with/ this program.  If not, see <https://www.gnu.org/licenses/>.
#ifndef GADDLE_BACKEND_GLOBAL_H
#define GADDLE_BACKEND_GLOBAL_H

#include "minimize.h"

#if defined(GADDLE_BACKEND_LIBRARY)
#  define GADDLE_BACKENDSHARED_EXPORT Q_DECL_EXPORT
#else
#  define GADDLE_BACKENDSHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // GADDLE_BACKEND_GLOBAL_H
