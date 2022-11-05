function castro2m2(; solver = nothing)
    m = JuMP.Model(solver)

    @variable(m, 0 <= x[1:41] <= 1E5)
    @variable(m, 0 <= obj <= 1E3)
    @objective(m, Min, obj)

    @NLconstraint(m, e31, x[28] * x[30] - x[16] == 0.0)
    @NLconstraint(m, e32, x[28] * x[31] - x[17] == 0.0)
    @NLconstraint(m, e33, x[29] * x[32] - x[18] == 0.0)
    @NLconstraint(m, e34, x[29] * x[33] - x[19] == 0.0)
    @NLconstraint(m, e35, x[28] * x[36] - x[22] == 0.0)
    @NLconstraint(m, e36, x[29] * x[37] - x[23] == 0.0)
    @NLconstraint(m, e37, x[14] * x[30] - x[1] == 0.0)
    @NLconstraint(m, e38, x[14] * x[31] - x[2] == 0.0)
    @NLconstraint(m, e39, x[15] * x[32] - x[3] == 0.0)
    @NLconstraint(m, e40, x[15] * x[33] - x[4] == 0.0)
    @NLconstraint(m, e41, x[14] * x[36] - x[7] == 0.0)
    @NLconstraint(m, e42, x[15] * x[37] - x[8] == 0.0)

    @constraint(m, e1, -x[14] - x[15] + obj == 0.0)
    @constraint(m, e2, -x[5] - x[9] - x[10] == -60.0)
    @constraint(m, e3, -x[6] - x[11] - x[12] == -20.0)
    @constraint(m, e4, -x[1] - x[3] - x[9] - x[11] + x[14] == 0.0)
    @constraint(m, e5, -x[2] - x[4] - x[10] - x[12] + x[15] == 0.0)
    @constraint(m, e6, -x[1] - x[2] - x[7] + x[14] == 0.0)
    @constraint(m, e7, -x[3] - x[4] - x[8] + x[15] == 0.0)
    @constraint(m, e8, -x[5] - x[6] - x[7] - x[8] + x[13] == 0.0)
    @constraint(m, e9, -x[20] - x[24] - x[25] == -24000.0)
    @constraint(m, e10, -x[21] - x[26] - x[27] == -16000.0)
    @constraint(m, e11, -x[24] + 24000 * x[38] == 0.0)
    @constraint(m, e12, -x[25] + 24000 * x[39] == 0.0)
    @constraint(m, e13, -x[26] + 16000 * x[40] == 0.0)
    @constraint(m, e14, -x[27] + 16000 * x[41] == 0.0)
    @constraint(m, e15, -x[20] + 24000 * x[34] == 0.0)
    @constraint(m, e16, -x[21] + 16000 * x[35] == 0.0)
    @constraint(m, e17, -x[9] + 60 * x[38] == 0.0)
    @constraint(m, e18, -x[10] + 60 * x[39] == 0.0)
    @constraint(m, e19, -x[11] + 20 * x[40] == 0.0)
    @constraint(m, e20, -x[12] + 20 * x[41] == 0.0)
    @constraint(m, e21, -x[5] + 60 * x[34] == 0.0)
    @constraint(m, e22, -x[6] + 20 * x[35] == 0.0)
    @constraint(m, e23, x[34] + x[38] + x[39] == 1.0)
    @constraint(m, e24, x[35] + x[40] + x[41] == 1.0)
    @constraint(m, e25, -200 * x[14] + x[16] + x[18] + x[24] + x[26] <= 0.0)
    @constraint(m, e26, -1000 * x[15] + x[17] + x[19] + x[25] + x[27] <= 0.0)
    @constraint(
        m,
        e27,
        0.01 * x[16] + 0.01 * x[18] + 0.01 * x[24] + 0.01 * x[26] - x[28] == 0.0
    )
    @constraint(
        m,
        e28,
        0.2 * x[17] + 0.2 * x[19] + 0.2 * x[25] + 0.2 * x[27] - x[29] == 0.0
    )
    @constraint(m, e29, -x[16] - x[17] - x[22] + x[28] == 0.0)
    @constraint(m, e30, -x[18] - x[19] - x[23] + x[29] == 0.0)
    @constraint(m, e43, x[30] + x[31] + x[36] == 1.0)
    @constraint(m, e44, x[32] + x[33] + x[37] == 1.0)
    @constraint(m, e45, -10 * x[13] + x[20] + x[21] + x[22] + x[23] <= 0.0)

    return m
end

function castro6m2(; solver = nothing)
    m = JuMP.Model(solver)

    @variable(m, 0 <= x[1:133] <= 1E6)
    @variable(m, obj)
    @variable(m, b[1:40])   #Additional variables to increase degree of freedom and Ipopt's convergence
    @objective(m, Min, obj)

    # Non-Linear Constraints
    @NLconstraint(m, e99, x[101] * x[110] - x[29] == 0.0)
    @NLconstraint(m, e100, x[102] * x[110] - x[30] == 0.0)
    @NLconstraint(m, e101, x[103] * x[110] - x[31] == 0.0)
    @NLconstraint(m, e102, x[101] * x[111] - x[32] == 0.0)
    @NLconstraint(m, e103, x[102] * x[111] - x[33] == 0.0)
    @NLconstraint(m, e104, x[103] * x[111] - x[34] == 0.0)
    @NLconstraint(m, e105, x[101] * x[112] - x[35] == 0.0)
    @NLconstraint(m, e106, x[102] * x[112] - x[36] == 0.0)
    @NLconstraint(m, e107, x[103] * x[112] - x[37] == 0.0)
    @NLconstraint(m, e108, x[104] * x[113] - x[38] == 0.0)
    @NLconstraint(m, e109, x[105] * x[113] - x[39] == 0.0)
    @NLconstraint(m, e110, x[106] * x[113] - x[40] == 0.0)
    @NLconstraint(m, e111, x[104] * x[114] - x[41] == 0.0)
    @NLconstraint(m, e112, x[105] * x[114] - x[42] == 0.0)
    @NLconstraint(m, e113, x[106] * x[114] - x[43] == 0.0)
    @NLconstraint(m, e114, x[104] * x[115] - x[44] == 0.0)
    @NLconstraint(m, e115, x[105] * x[115] - x[45] == 0.0)
    @NLconstraint(m, e116, x[106] * x[115] - x[46] == 0.0)
    @NLconstraint(m, e117, x[107] * x[116] - x[47] == 0.0)
    @NLconstraint(m, e118, x[108] * x[116] - x[48] == 0.0)
    @NLconstraint(m, e119, x[109] * x[116] - x[49] == 0.0)
    @NLconstraint(m, e120, x[107] * x[117] - x[50] == 0.0)
    @NLconstraint(m, e121, x[108] * x[117] - x[51] == 0.0)
    @NLconstraint(m, e122, x[109] * x[117] - x[52] == 0.0)
    @NLconstraint(m, e123, x[107] * x[118] - x[53] == 0.0)
    @NLconstraint(m, e124, x[108] * x[118] - x[54] == 0.0)
    @NLconstraint(m, e125, x[109] * x[118] - x[55] == 0.0)
    @NLconstraint(m, e126, x[101] * x[122] - x[65] == 0.0)
    @NLconstraint(m, e127, x[102] * x[122] - x[66] == 0.0)
    @NLconstraint(m, e128, x[103] * x[122] - x[67] == 0.0)
    @NLconstraint(m, e129, x[104] * x[123] - x[68] == 0.0)
    @NLconstraint(m, e130, x[105] * x[123] - x[69] == 0.0)
    @NLconstraint(m, e131, x[106] * x[123] - x[70] == 0.0)
    @NLconstraint(m, e132, x[107] * x[124] - x[71] == 0.0)
    @NLconstraint(m, e133, x[108] * x[124] - x[72] == 0.0)
    @NLconstraint(m, e134, x[109] * x[124] - x[73] == 0.0)
    @NLconstraint(m, e135, x[26] * x[110] - x[1] == 0.0)
    @NLconstraint(m, e136, x[26] * x[111] - x[2] == 0.0)
    @NLconstraint(m, e137, x[26] * x[112] - x[3] == 0.0)
    @NLconstraint(m, e138, x[27] * x[113] - x[4] == 0.0)
    @NLconstraint(m, e139, x[27] * x[114] - x[5] == 0.0)
    @NLconstraint(m, e140, x[27] * x[115] - x[6] == 0.0)
    @NLconstraint(m, e141, x[28] * x[116] - x[7] == 0.0)
    @NLconstraint(m, e142, x[28] * x[117] - x[8] == 0.0)
    @NLconstraint(m, e143, x[28] * x[118] - x[9] == 0.0)
    @NLconstraint(m, e144, x[26] * x[122] - x[13] == 0.0)
    @NLconstraint(m, e145, x[27] * x[123] - x[14] == 0.0)
    @NLconstraint(m, e146, x[28] * x[124] - x[15] == 0.0)

    @constraint(m, e1, -x[26] - x[27] - x[28] + obj == 0.0)
    @constraint(m, e2, -x[10] - x[16] - x[17] - x[18] == -13.1)
    @constraint(m, e3, -x[11] - x[19] - x[20] - x[21] == -32.7)
    @constraint(m, e4, -x[12] - x[22] - x[23] - x[24] == -56.5)
    @constraint(m, e5, -x[1] - x[4] - x[7] - x[16] - x[19] - x[22] + x[26] == 0.0)
    @constraint(m, e6, -x[2] - x[5] - x[8] - x[17] - x[20] - x[23] + x[27] == 0.0)
    @constraint(m, e7, -x[3] - x[6] - x[9] - x[18] - x[21] - x[24] + x[28] == 0.0)
    @constraint(m, e8, -x[1] - x[2] - x[3] - x[13] + x[26] == 0.0)
    @constraint(m, e9, -x[4] - x[5] - x[6] - x[14] + x[27] == 0.0)
    @constraint(m, e10, -x[7] - x[8] - x[9] - x[15] + x[28] == 0.0)
    @constraint(m, e11, -x[10] - x[11] - x[12] - x[13] - x[14] - x[15] + x[25] == 0.0)
    @constraint(m, e12, -x[56] - x[74] - x[77] - x[80] == -131.0)
    @constraint(m, e13, -x[57] - x[75] - x[78] - x[81] == -5109.0)
    @constraint(m, e14, -x[58] - x[76] - x[79] - x[82] == -3275.0)
    @constraint(m, e15, -x[59] - x[83] - x[86] - x[89] == -3597.0)
    @constraint(m, e16, -x[60] - x[84] - x[87] - x[90] == -548706.0)
    @constraint(m, e17, -x[61] - x[85] - x[88] - x[91] == -13080.0)
    @constraint(m, e18, -x[62] - x[92] - x[95] - x[98] == -5650.0)
    @constraint(m, e19, -x[63] - x[93] - x[96] - x[99] == -1412.5)
    @constraint(m, e20, -x[64] - x[94] - x[97] - x[100] == -19775.0)
    @constraint(m, e21, -x[74] + 131 * x[125] == 0.0)
    @constraint(m, e22, -x[75] + 5109 * x[125] == 0.0)
    @constraint(m, e23, -x[76] + 3275 * x[125] == 0.0)
    @constraint(m, e24, -x[77] + 131 * x[126] == 0.0)
    @constraint(m, e25, -x[78] + 5109 * x[126] == 0.0)
    @constraint(m, e26, -x[79] + 3275 * x[126] == 0.0)
    @constraint(m, e27, -x[80] + 131 * x[127] == 0.0)
    @constraint(m, e28, -x[81] + 5109 * x[127] == 0.0)
    @constraint(m, e29, -x[82] + 3275 * x[127] == 0.0)
    @constraint(m, e30, -x[83] + 3597 * x[128] == 0.0)
    @constraint(m, e31, -x[84] + 548706 * x[128] == 0.0)
    @constraint(m, e32, -x[85] + 13080 * x[128] == 0.0)
    @constraint(m, e33, -x[86] + 3597 * x[129] == 0.0)
    @constraint(m, e34, -x[87] + 548706 * x[129] == 0.0)
    @constraint(m, e35, -x[88] + 13080 * x[129] == 0.0)
    @constraint(m, e36, -x[89] + 3597 * x[130] == 0.0)
    @constraint(m, e37, -x[90] + 548706 * x[130] == 0.0)
    @constraint(m, e38, -x[91] + 13080 * x[130] == 0.0)
    @constraint(m, e39, -x[92] + 5650 * x[131] == 0.0)
    @constraint(m, e40, -x[93] + 1412.5 * x[131] == 0.0)
    @constraint(m, e41, -x[94] + 19775 * x[131] == 0.0)
    @constraint(m, e42, -x[95] + 5650 * x[132] == 0.0)
    @constraint(m, e43, -x[96] + 1412.5 * x[132] == 0.0)
    @constraint(m, e44, -x[97] + 19775 * x[132] == 0.0)
    @constraint(m, e45, -x[98] + 5650 * x[133] == 0.0)
    @constraint(m, e46, -x[99] + 1412.5 * x[133] == 0.0)
    @constraint(m, e47, -x[100] + 19775 * x[133] == 0.0)
    @constraint(m, e48, -x[56] + 131 * x[119] == 0.0)
    @constraint(m, e49, -x[57] + 5109 * x[119] == 0.0)
    @constraint(m, e50, -x[58] + 3275 * x[119] == 0.0)
    @constraint(m, e51, -x[59] + 3597 * x[120] == 0.0)
    @constraint(m, e52, -x[60] + 548706 * x[120] == 0.0)
    @constraint(m, e53, -x[61] + 13080 * x[120] == 0.0)
    @constraint(m, e54, -x[62] + 5650 * x[121] == 0.0)
    @constraint(m, e55, -x[63] + 1412.5 * x[121] == 0.0)
    @constraint(m, e56, -x[64] + 19775 * x[121] == 0.0)
    @constraint(m, e57, -x[16] + 13.1 * x[125] == 0.0)
    @constraint(m, e58, -x[17] + 13.1 * x[126] == 0.0)
    @constraint(m, e59, -x[18] + 13.1 * x[127] == 0.0)
    @constraint(m, e60, -x[19] + 32.7 * x[128] == 0.0)
    @constraint(m, e61, -x[20] + 32.7 * x[129] == 0.0)
    @constraint(m, e62, -x[21] + 32.7 * x[130] == 0.0)
    @constraint(m, e63, -x[22] + 56.5 * x[131] == 0.0)
    @constraint(m, e64, -x[23] + 56.5 * x[132] == 0.0)
    @constraint(m, e65, -x[24] + 56.5 * x[133] == 0.0)
    @constraint(m, e66, -x[10] + 13.1 * x[119] == 0.0)
    @constraint(m, e67, -x[11] + 32.7 * x[120] == 0.0)
    @constraint(m, e68, -x[12] + 56.5 * x[121] == 0.0)
    @constraint(m, e69, x[119] + x[125] + x[126] + x[127] == 1.0)
    @constraint(m, e70, x[120] + x[128] + x[129] + x[130] == 1.0)
    @constraint(m, e71, x[121] + x[131] + x[132] + x[133] == 1.0)
    @constraint(
        m,
        e72,
        -110 * x[26] + x[29] + x[38] + x[47] + x[74] + x[83] + x[92] <= 0.0
    )
    @constraint(
        m,
        e73,
        -16780 * x[26] + x[30] + x[39] + x[48] + x[75] + x[84] + x[93] <= 0.0
    )
    @constraint(
        m,
        e74,
        -400 * x[26] + x[31] + x[40] + x[49] + x[76] + x[85] + x[94] <= 0.0
    )
    @constraint(
        m,
        e75,
        -110 * x[27] + x[32] + x[41] + x[50] + x[77] + x[86] + x[95] <= 0.0
    )
    @constraint(
        m,
        e76,
        -16780 * x[27] + x[33] + x[42] + x[51] + x[78] + x[87] + x[96] <= 0.0
    )
    @constraint(
        m,
        e77,
        -400 * x[27] + x[34] + x[43] + x[52] + x[79] + x[88] + x[97] <= 0.0
    )
    @constraint(
        m,
        e78,
        -110 * x[28] + x[35] + x[44] + x[53] + x[80] + x[89] + x[98] <= 0.0
    )
    @constraint(
        m,
        e79,
        -16780 * x[28] + x[36] + x[45] + x[54] + x[81] + x[90] + x[99] <= 0.0
    )
    @constraint(
        m,
        e80,
        -400 * x[28] + x[37] + x[46] + x[55] + x[82] + x[91] + x[100] <= 0.0
    )
    @constraint(m, e81, x[29] + x[38] + x[47] + x[74] + x[83] + x[92] - x[101] == 0.0)
    @constraint(
        m,
        e82,
        0.001 * x[30] +
        0.001 * x[39] +
        0.001 * x[48] +
        0.001 * x[75] +
        0.001 * x[84] +
        0.001 * x[93] - x[102] == 0.0
    )
    @constraint(m, e83, x[31] + x[40] + x[49] + x[76] + x[85] + x[94] - x[103] == 0.0)
    @constraint(
        m,
        e84,
        0.3 * x[32] +
        0.3 * x[41] +
        0.3 * x[50] +
        0.3 * x[77] +
        0.3 * x[86] +
        0.3 * x[95] - x[104] == 0.0
    )
    @constraint(
        m,
        e85,
        0.1 * x[33] +
        0.1 * x[42] +
        0.1 * x[51] +
        0.1 * x[78] +
        0.1 * x[87] +
        0.1 * x[96] - x[105] == 0.0
    )
    @constraint(
        m,
        e86,
        0.02 * x[34] +
        0.02 * x[43] +
        0.02 * x[52] +
        0.02 * x[79] +
        0.02 * x[88] +
        0.02 * x[97] - x[106] == 0.0
    )
    @constraint(
        m,
        e87,
        0.3 * x[35] +
        0.3 * x[44] +
        0.3 * x[53] +
        0.3 * x[80] +
        0.3 * x[89] +
        0.3 * x[98] - x[107] == 0.0
    )
    @constraint(m, e88, x[36] + x[45] + x[54] + x[81] + x[90] + x[99] - x[108] == 0.0)
    @constraint(
        m,
        e89,
        0.5 * x[37] +
        0.5 * x[46] +
        0.5 * x[55] +
        0.5 * x[82] +
        0.5 * x[91] +
        0.5 * x[100] - x[109] == 0.0
    )
    @constraint(m, e90, -x[29] - x[32] - x[35] - x[65] + x[101] == 0.0)
    @constraint(m, e91, -x[30] - x[33] - x[36] - x[66] + x[102] == 0.0)
    @constraint(m, e92, -x[31] - x[34] - x[37] - x[67] + x[103] == 0.0)
    @constraint(m, e93, -x[38] - x[41] - x[44] - x[68] + x[104] == 0.0)
    @constraint(m, e94, -x[39] - x[42] - x[45] - x[69] + x[105] == 0.0)
    @constraint(m, e95, -x[40] - x[43] - x[46] - x[70] + x[106] == 0.0)
    @constraint(m, e96, -x[47] - x[50] - x[53] - x[71] + x[107] == 0.0)
    @constraint(m, e97, -x[48] - x[51] - x[54] - x[72] + x[108] == 0.0)
    @constraint(m, e98, -x[49] - x[52] - x[55] - x[73] + x[109] == 0.0)
    @constraint(m, e147, x[110] + x[111] + x[112] + x[122] == 1.0)
    @constraint(m, e148, x[113] + x[114] + x[115] + x[123] == 1.0)
    @constraint(m, e149, x[116] + x[117] + x[118] + x[124] == 1.0)
    @constraint(
        m,
        e150,
        -20 * x[25] + x[56] + x[59] + x[62] + x[65] + x[68] + x[71] <= 0.0
    )
    @constraint(
        m,
        e151,
        -5 * x[25] + x[57] + x[60] + x[63] + x[66] + x[69] + x[72] <= 0.0
    )
    @constraint(
        m,
        e152,
        -100 * x[25] + x[58] + x[61] + x[64] + x[67] + x[70] + x[73] <= 0.0
    )

    return m
end

function castro4m2(; solver = nothing)
    m = JuMP.Model(solver)

    # ----- Variables ----- #
    @variable(m, 0 <= x[1:55] <= 1E6)
    @variable(m, obj)
    @variable(m, b[1:20]) #Additional variables to increase degree of freedom and Ipopt's convergence
    @objective(m, Min, obj)

    @NLconstraint(m, e45, x[40] * x[44] - x[16] == 0.0)
    @NLconstraint(m, e46, x[41] * x[44] - x[17] == 0.0)
    @NLconstraint(m, e47, x[40] * x[45] - x[18] == 0.0)
    @NLconstraint(m, e48, x[41] * x[45] - x[19] == 0.0)
    @NLconstraint(m, e49, x[42] * x[46] - x[20] == 0.0)
    @NLconstraint(m, e50, x[43] * x[46] - x[21] == 0.0)
    @NLconstraint(m, e51, x[42] * x[47] - x[22] == 0.0)
    @NLconstraint(m, e52, x[43] * x[47] - x[23] == 0.0)
    @NLconstraint(m, e53, x[40] * x[50] - x[28] == 0.0)
    @NLconstraint(m, e54, x[41] * x[50] - x[29] == 0.0)
    @NLconstraint(m, e55, x[42] * x[51] - x[30] == 0.0)
    @NLconstraint(m, e56, x[43] * x[51] - x[31] == 0.0)
    @NLconstraint(m, e57, x[14] * x[44] - x[1] == 0.0)
    @NLconstraint(m, e58, x[14] * x[45] - x[2] == 0.0)
    @NLconstraint(m, e59, x[15] * x[46] - x[3] == 0.0)
    @NLconstraint(m, e60, x[15] * x[47] - x[4] == 0.0)
    @NLconstraint(m, e61, x[14] * x[50] - x[7] == 0.0)
    @NLconstraint(m, e62, x[15] * x[51] - x[8] == 0.0)

    @constraint(m, e1, -x[14] - x[15] + obj == 0.0)
    @constraint(m, e2, -x[5] - x[9] - x[10] == -40.0)
    @constraint(m, e3, -x[6] - x[11] - x[12] == -40.0)
    @constraint(m, e4, -x[1] - x[3] - x[9] - x[11] + x[14] == 0.0)
    @constraint(m, e5, -x[2] - x[4] - x[10] - x[12] + x[15] == 0.0)
    @constraint(m, e6, -x[1] - x[2] - x[7] + x[14] == 0.0)
    @constraint(m, e7, -x[3] - x[4] - x[8] + x[15] == 0.0)
    @constraint(m, e8, -x[5] - x[6] - x[7] - x[8] + x[13] == 0.0)
    @constraint(m, e9, -x[24] - x[32] - x[34] == -4000.0)
    @constraint(m, e10, -x[25] - x[33] - x[35] == -800.0)
    @constraint(m, e11, -x[26] - x[36] - x[38] == -600.0)
    @constraint(m, e12, -x[27] - x[37] - x[39] == -8000.0)
    @constraint(m, e13, -x[32] + 4000 * x[52] == 0.0)
    @constraint(m, e14, -x[33] + 800 * x[52] == 0.0)
    @constraint(m, e15, -x[34] + 4000 * x[53] == 0.0)
    @constraint(m, e16, -x[35] + 800 * x[53] == 0.0)
    @constraint(m, e17, -x[36] + 600 * x[54] == 0.0)
    @constraint(m, e18, -x[37] + 8000 * x[54] == 0.0)
    @constraint(m, e19, -x[38] + 600 * x[55] == 0.0)
    @constraint(m, e20, -x[39] + 8000 * x[55] == 0.0)
    @constraint(m, e21, -x[24] + 4000 * x[48] == 0.0)
    @constraint(m, e22, -x[25] + 800 * x[48] == 0.0)
    @constraint(m, e23, -x[26] + 600 * x[49] == 0.0)
    @constraint(m, e24, -x[27] + 8000 * x[49] == 0.0)
    @constraint(m, e25, -x[9] + 40 * x[52] == 0.0)
    @constraint(m, e26, -x[10] + 40 * x[53] == 0.0)
    @constraint(m, e27, -x[11] + 40 * x[54] == 0.0)
    @constraint(m, e28, -x[12] + 40 * x[55] == 0.0)
    @constraint(m, e29, -x[5] + 40 * x[48] == 0.0)
    @constraint(m, e30, -x[6] + 40 * x[49] == 0.0)
    @constraint(m, e31, x[48] + x[52] + x[53] == 1.0)
    @constraint(m, e32, x[49] + x[54] + x[55] == 1.0)
    @constraint(m, e33, -200 * x[14] + x[16] + x[20] + x[32] + x[36] <= 0.0)
    @constraint(m, e34, -200 * x[14] + x[17] + x[21] + x[33] + x[37] <= 0.0)
    @constraint(m, e35, -200 * x[15] + x[18] + x[22] + x[34] + x[38] <= 0.0)
    @constraint(m, e36, -200 * x[15] + x[19] + x[23] + x[35] + x[39] <= 0.0)
    @constraint(
        m,
        e37,
        0.05 * x[16] + 0.05 * x[20] + 0.05 * x[32] + 0.05 * x[36] - x[40] == 0.0
    )
    @constraint(m, e38, x[17] + x[21] + x[33] + x[37] - x[41] == 0.0)
    @constraint(m, e39, x[18] + x[22] + x[34] + x[38] - x[42] == 0.0)
    @constraint(
        m,
        e40,
        0.024 * x[19] + 0.024 * x[23] + 0.024 * x[35] + 0.024 * x[39] - x[43] == 0.0
    )
    @constraint(m, e41, -x[16] - x[18] - x[28] + x[40] == 0.0)
    @constraint(m, e42, -x[17] - x[19] - x[29] + x[41] == 0.0)
    @constraint(m, e43, -x[20] - x[22] - x[30] + x[42] == 0.0)
    @constraint(m, e44, -x[21] - x[23] - x[31] + x[43] == 0.0)
    @constraint(m, e63, x[44] + x[45] + x[50] == 1.0)
    @constraint(m, e64, x[46] + x[47] + x[51] == 1.0)
    @constraint(m, e65, -10 * x[13] + x[24] + x[26] + x[28] + x[30] <= 0.0)
    @constraint(m, e66, -10 * x[13] + x[25] + x[27] + x[29] + x[31] <= 0.0)

    return m
end
