//
//  Model.swift
//  calcGUI2
//
//  Created by (s) Mickael Hinman on 18/02/2016.
//  Copyright (c) 2016 (s) Mickael Hinman. All rights reserved.
//

import Foundation
import Accelerate

func addCC(constant: Float, constant2: Float) -> Float{                                         //Constant add Constant
    let total = constant + constant2
    return total
}

func addCCZ(constant: [Float], constant2: [Float]) -> [Float]{                                  //Constant add Constant Complex
    let total: [Float] = [constant[0] + constant2[0], constant[1] + constant2[1]]
    return total
}

func subtractCC(constant: Float, constant2: Float) -> Float{                                    //Constant subtract Constant
    let total = constant - constant2
    return total
}

func subtractCCZ(constant: [Float], constant2: [Float]) -> [Float]{                             //Constant subtract Constant Complex
    let total: [Float] = [constant[0] - constant2[0], constant[1] - constant2[1]]
    return total
}

func multiplyCC(constant: Float, constant2: Float) -> Float{                                    //Constant multiply Constant
    let total = constant * constant2
    return total
}

func multiplyCCZ(constant: [Float], constant2: [Float]) -> [Float]{                             //Constant multiply Constant Complex
    let a: Float = constant[0]
    let b: Float = constant[1]
    let c: Float = constant2[0]
    let d: Float = constant2[1]
    
    let realConstant = a * c - b * d
    let imaginaryConstant = a * d + b * c
    
    let total: [Float] = [realConstant, imaginaryConstant]
    return total
}

func divideCC(constant: Float, constant2: Float) -> Float{                                      //Constant divide Constant
    var total: Float = 0
    if(constant2 != 0){
        total = constant / constant2
    }
    else{
        print("error divideCC: dividing by 0")
    }
    return total
}

func divideCCZ(constant: [Float], constant2: [Float]) -> [Float]{                               //Constant divide Constant Complex
    let a: Float = constant[0]
    let b: Float = constant[1]
    let c: Float = constant2[0]
    let d: Float = constant2[1]
    
    let realConstant = (a * c + b * d) / (c * c + d * d)
    var imaginaryConstant: Float = 0.0
    if(!(c == 0 && d == 0)){
        imaginaryConstant = (b * c - a * d) / (c * c + d * d)
    }
    else{
        print("error divideCCZ: dividing by 0")
    }
    
    let total: [Float] = [realConstant, imaginaryConstant]
    return total
}

func addVC(vector: [Float], var constant: Float) -> [Float]{                                    //Vector add Constant
    var total: [Float] = vector
    vDSP_vsadd(vector, 1, &constant, &total, 1, UInt(vector.count))
    return total
}

func addVCZ(vector: [Float], constant: [Float]) -> [Float]{                                     //Vector add Constant Complex
    var total = [Float](count: vector.count, repeatedValue: 0.0)
    
    for var index = 0; index < vector.count; index += 2{
        let temp = addCCZ([vector[index], vector[index + 1]], [constant[0], constant[1]])
        total[index] = temp[0]
        total[index + 1] = temp[1]
    }
    return total
}

func subtractVC(vector: [Float], constant: Float) -> [Float]{                                   //Vector subtract Constant
    var negatedConstant = -constant
    var total: [Float] = vector
    vDSP_vsadd(vector, 1, &negatedConstant, &total, 1, UInt(vector.count))
    return total
}

func subtractVCZ(vector: [Float], constant: [Float]) -> [Float]{                                //Vector subtract Constant Complex
    var total = [Float](count: vector.count, repeatedValue: 0.0)
    
    for var index = 0; index < vector.count; index += 2{
        let temp = subtractCCZ([vector[index], vector[index + 1]], [constant[0], constant[1]])
        total[index] = temp[0]
        total[index + 1] = temp[1]
    }
    return total
}

func multiplyVC(vector: [Float], var constant: Float) -> [Float]{                               //Vector multiply Constant
    var total: [Float] = vector
    vDSP_vsmul(vector, 1, &constant, &total, 1, UInt(vector.count))
    return total
}

func multiplyVCZ(vector: [Float], var constant: [Float]) -> [Float]{                            //Vector multiply Constant Complex
    var total = [Float](count: vector.count, repeatedValue: 0.0)
    
    for var index = 0; index < vector.count; index += 2{
        let temp = multiplyCCZ([vector[index], vector[index + 1]], [constant[0], constant[1]])
        total[index] = temp[0]
        total[index + 1] = temp[1]
    }
    return total
}

func divideVC(vector: [Float], constant: Float) -> [Float]{                                     //Vector divide Constant
    var dividedConstant = 1 / constant
    var total: [Float] = vector
    vDSP_vsmul(vector, 1, &dividedConstant, &total, 1, UInt(vector.count))
    return total
}

func divideVCZ(vector: [Float], constant: [Float]) -> [Float]{                                  //Vector divide Constant Complex
    var total = [Float](count: vector.count, repeatedValue: 0.0)
    
    for var index = 0; index < vector.count; index += 2{
        let temp = divideCCZ([vector[index], vector[index + 1]], [constant[0], constant[1]])
        total[index] = temp[0]
        total[index + 1] = temp[1]
    }
    return total
}

func addVV(vector: [Float], vector2: [Float]) -> [Float]{                                       //Vector add Vector
    var total: [Float] = vector
    vDSP_vadd(vector, 1, vector2, 1, &total, 1, UInt(vector.count))
    return total
}

func addVVZ(vector: [Float], vector2: [Float]) -> [Float]{                                      //Vector add Vector Complex
    var total = [Float](count: vector.count, repeatedValue: 0.0)
    
    for var index = 0; index < vector.count; index += 2{
        let temp = addCCZ([vector[index], vector[index + 1]], [vector2[index], vector2[index + 1]])
        total[index] = temp[0]
        total[index + 1] = temp[1]
    }
    return total
}

func subtractVV(vector: [Float], vector2: [Float]) -> [Float]{                                  //Vector subtract Vector
    var vector2Negated: [Float] = vector2
    let vector2NegatedLength = vector2Negated.count
    for index in 0..<vector2NegatedLength{
        vector2Negated[index] = -vector2Negated[index]
    }
    var total: [Float] = vector
    vDSP_vadd(vector, 1, vector2Negated, 1, &total, 1, UInt(vector.count))
    return total
}

func subtractVVZ(vector: [Float], vector2: [Float]) -> [Float]{                                 //Vector subtract Vector Complex
    var total = [Float](count: vector.count, repeatedValue: 0.0)
    
    for var index = 0; index < vector.count; index += 2{
        let temp = subtractCCZ([vector[index], vector[index + 1]], [vector2[index], vector2[index + 1]])
        total[index] = temp[0]
        total[index + 1] = temp[1]
    }
    return total
}

func addMC(matrix: [Float], constant: Float) -> [Float]{                                        //Matrix add Constant
    return addMM(matrix, multiplyMC(identityMatrix(Int(sqrt(Float(matrix.count)))).matrix, constant))
}

func addMCZ(matrix: [Float], constant: [Float]) -> [Float]{                                     //Matrix add Constant Complex
    return addVCZ(matrix, constant)
}

func subtractMC(matrix: [Float], constant: Float) -> [Float]{                                   //Matrix subtract Constant
    return subtractVC(matrix, constant)
}

func subtractMCZ(matrix: [Float], constant: [Float]) -> [Float]{                                //Matrix subtract Constant Complex
    return subtractVCZ(matrix, constant)
}

func multiplyMC(matrix: [Float], constant: Float) -> [Float]{                                   //Matrix multiply Constant
    return multiplyVC(matrix, constant)
}

func multiplyMCZ(matrix: [Float], constant: [Float]) -> [Float]{                                //Matrix multiply Constant Complex
    return multiplyVCZ(matrix, constant)
}

func divideMC(matrix: [Float], constant: Float) -> [Float]{                                     //Matrix divide Constant
    return divideVC(matrix, constant)
}

func divideMCZ(matrix: [Float], constant: [Float]) -> [Float]{                                  //Matrix divide Constant Complex
    return divideVCZ(matrix, constant)
}

func multiplyMV(matrix: [Float], vector: [Float]) -> [Float]{                                   //Matrix multiply Vector
    let matrixSize = sqrt(Double(matrix.count))
    var total: [Float] = [Float](count: Int(matrixSize), repeatedValue: 0.0)
    cblas_sgemv(CblasRowMajor, CblasNoTrans, Int32(matrixSize), Int32(matrixSize), 1.0, matrix, Int32(matrixSize), vector, 1, 0, &total, 1)
    return total
}

func multiplyMVZ(matrix: [Float], vector: [Float]) -> [Float]{                                  //Matrix multiply Vector Complex
    let matrixSize = sqrt(Double(matrix.count / 2))
    let alpha: [Float] = [1.0, 0.0]
    let beta: [Float] = [0.0, 0.0]
    var total: [Float] = [Float](count: vector.count, repeatedValue: 0.0)
    cblas_cgemv(CblasRowMajor, CblasNoTrans, Int32(matrixSize), Int32(matrixSize), alpha, matrix, Int32(matrixSize), vector, 1, beta, &total, 1)
    return total
}

func addMM(matrix: [Float], matrix2: [Float]) -> [Float]{                                       //Matrix add Matrix
    return addVV(matrix, matrix2)
}

func addMMZ(matrix: [Float], matrix2: [Float]) -> [Float]{                                      //Matrix add Matrix Complex
    return addVVZ(matrix, matrix2)
}

func subtractMM(matrix: [Float], matrix2: [Float]) -> [Float]{                                  //Matrix subtract Matrix
    return subtractVV(matrix, matrix2)
}

func subtractMMZ(matrix: [Float], matrix2: [Float]) -> [Float]{                                 //Matrix subtract Matrix Complex
    return subtractVVZ(matrix, matrix2)
}

func multiplyMM(matrix: [Float], matrix2: [Float]) -> [Float]{                                  //Matrix multiply Matrix
    let matrixSize = sqrt(Double(matrix.count))
    var total: [Float] = [Float](count: matrix.count, repeatedValue: 0.0)
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Int32(matrixSize), Int32(matrixSize), Int32(matrixSize), 1.0, matrix, Int32(matrixSize), matrix2, Int32(matrixSize), 0, &total, Int32(matrixSize))
    return total
}

func multiplyMMZ(matrix: [Float], matrix2: [Float]) -> [Float]{                                 //Matrix multiply Matrix Complex
    let matrixSize = sqrt(Double(matrix.count / 2))
    let alpha: [Float] = [1.0, 0.0]
    let beta: [Float] = [0.0, 0.0]
    var total: [Float] = [Float](count: matrix.count, repeatedValue: 0.0)
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Int32(matrixSize), Int32(matrixSize), Int32(matrixSize), alpha, matrix, Int32(matrixSize), matrix2, Int32(matrixSize), beta, &total, Int32(matrixSize))
    return total
}

func negateC(constant: Float) -> Float{                                                           //Constant negate
    return -constant
}

func negateCZ(constant: [Float]) -> [Float]{                                                      //Constant negate Complex
    return [-constant[0], -constant[1]]
}

func negateV(vector: [Float]) -> [Float]{                                                         //Vector negate
    var total: [Float] = [Float](count: vector.count, repeatedValue: 0.0)
    for(var index = 0; index < vector.count; index += 1){
        total[index] = -vector[index]
    }
    return total
}

func negateVZ(vector: [Float]) -> [Float]{                                                        //Vector negate Complex
    return negateV(vector)
}

func negateM(matrix: [Float]) -> [Float]{                                                         //Matrix negate
    return negateV(matrix)
}

func negateMZ(matrix: [Float]) -> [Float]{                                                        //Matrix negate Complex
    return negateV(matrix)
}

func transposeM(matrix: [Float]) -> [Float]{
    var total: [Float] = [Float](count: matrix.count, repeatedValue: 0.0)
    let matrixSize = Int(sqrt(Double(matrix.count)))
    for index in 0..<matrixSize{
        for number in 0..<matrixSize{
            total[index * matrixSize + number] = matrix[index + number * matrixSize]
        }
    }
    return total
}

func transposeMZ(matrix: [Float]) -> [Float]{
    var total: [Float] = [Float](count: matrix.count, repeatedValue: 0.0)
    let matrixSize = Int(sqrt(Double(matrix.count / 2)))
    for index in 0..<matrixSize{
        for number in 0..<matrixSize{
            total[index * matrixSize + number] = matrix[2 * index + 2 * number * matrixSize]
            total[index * matrixSize + number] = matrix[2 * index + 1 + 2 * number * matrixSize]
        }
    }
    return total
}

func complexConjugateCZ(constant: [Float]) -> [Float]{
    return [constant[0], -constant[1]]
}

func complexConjugateVZ(vector: [Float]) -> [Float]{
    var total: [Float] = [Float](count: vector.count, repeatedValue: 0.0)
    for(var index = 0; index < vector.count; index += 2){
        total[index] = vector[index]
        total[index + 1] = -vector[index + 1]
    }
    return total
}

func complexConjugateMZ(matrix: [Float]) -> [Float]{
    return complexConjugateVZ(matrix)
}

func daggerM(matrix: [Float]) -> [Float]{
    return transposeM(matrix)
}

func daggerMZ(matrix: [Float]) -> [Float]{
    return complexConjugateMZ(transposeMZ(matrix))
}

func squaredC(constant: Float) -> Float{
    return multiplyCC(constant, constant)
}

func squaredCZ(constant: [Float]) -> [Float]{
    return multiplyCCZ(constant, constant)
}

func squaredM(matrix: [Float]) -> [Float]{
    return multiplyMM(matrix, matrix)
}

func squaredMZ(matrix: [Float]) -> [Float]{
    return multiplyMMZ(matrix, matrix)
}

func cubedC(constant: Float) -> Float{
    return multiplyCC(multiplyCC(constant, constant), constant)
}

func cubedCZ(constant: [Float]) -> [Float]{
    return multiplyCCZ(multiplyCCZ(constant, constant), constant)
}

func cubedM(matrix: [Float]) -> [Float]{
    return multiplyMM(multiplyMM(matrix, matrix), matrix)
}

func cubedMZ(matrix: [Float]) -> [Float]{
    return multiplyMMZ(multiplyMMZ(matrix, matrix), matrix)
}

func squareRootC(constant: Float) -> Float{
    return Float(sqrt(Double(constant)))
}

class EquationElement {
    enum types{
        case VECTOR
        case VECTORCOMPLEX
        case MATRIX
        case MATRIXCOMPLEX
        case BRACKET
        case OPERATION
        case CONSTANT
        case CONSTANTCOMPLEX
        case NOTSET
    }
    enum bracketTypes{
        case NONE
        case LEFT
        case RIGHT
    }
    enum operatorTypes{
        case NONE
        case ADD
        case SUBTRACT
        case MULTIPLY
        case DIVIDE
        case TRANSPOSE
        case COMPLEXCONJUGATE
        case DAGGER
        case NEGATE
        case POWER
        case SQUARED
        case CUBED
        case SQUAREROOT
        case FACTORIAL
        case COS
        case SIN
        case EXP
        case COSQ
    }
    
    var type: types = .NOTSET
    var vector: [Float] = []
    var vectorComplex: [Float] = []
    var matrix: [Float] = []
    var matrixComplex: [Float] = []
    var braket: bracketTypes = .NONE
    var operation: operatorTypes = .NONE
    var constant: Float = 0
    var constantComplex: [Float] = []
    
    init(_ typeVar: types, _ arrayVar: [Float]){
        if(typeVar == .VECTOR){
            vector = arrayVar
            type = typeVar
        }
        else if(typeVar == .VECTORCOMPLEX){
            vectorComplex = arrayVar
            type = typeVar
        }
        else if(typeVar == .MATRIX){
            matrix = arrayVar
            type = typeVar
        }
        else if(typeVar == .MATRIXCOMPLEX){
            matrixComplex = arrayVar
            type = typeVar
        }
        else if(typeVar == .CONSTANTCOMPLEX){
            constantComplex = arrayVar
            type = typeVar
        }
    }
    init(_ bracketVar: bracketTypes){
        braket = bracketVar
        type = .BRACKET
    }
    init(_ operationVar: operatorTypes){
        operation = operationVar
        type = .OPERATION
    }
    init(_ constantVar: Float){
        constant = constantVar
        type = .CONSTANT
    }
    init(){
    }
    
}

func identityMatrix(matrixSize: Int) -> EquationElement{
    let result: EquationElement = EquationElement()
    result.type = .MATRIX
    result.matrix = [Float](count: matrixSize * matrixSize, repeatedValue: 0.0)
    
    for var row = 0; row < matrixSize; row += 1{
        for var column = 0; column < matrixSize; column += 1{
            if(row == column){
                result.matrix[row * matrixSize + column] = 1.0
            }
        }
    }
    return result
}

func absoluteValue(value: Float) -> Float{
    var result = value
    if(value < 0){
        result *= -1
    }
    return result
}

func maximumDifference(firstElement: EquationElement, secondElement: EquationElement) -> Float{
    var result: Float = 0
    let tempElement = implimentElementArraySimplification([firstElement, EquationElement(.SUBTRACT), secondElement])
    
    if(tempElement.type == .CONSTANT){
        let tempResult = absoluteValue(tempElement.constant)
        if(tempResult > result){
            result = tempResult
        }
    }
    else if(tempElement.type == .CONSTANTCOMPLEX){
        for(var index = 0; index < tempElement.constantComplex.count; index += 1){
            var tempResult = tempElement.constantComplex[index]
            tempResult = absoluteValue(tempResult)
            if(tempResult > result){
                result = tempResult
            }
        }
    }
    else if(tempElement.type == .MATRIX){
        for(var index = 0; index < tempElement.matrix.count; index += 1){
            var tempResult = tempElement.matrix[index]
            tempResult = absoluteValue(tempResult)
            if(tempResult > result){
                result = tempResult
            }
        }
    }
    else if(tempElement.type == .MATRIXCOMPLEX){
        for(var index = 0; index < tempElement.matrixComplex.count; index += 1){
            var tempResult = tempElement.matrixComplex[index]
            tempResult = absoluteValue(tempResult)
            if(tempResult > result){
                result = tempResult
            }
        }
    }
    return result
}

func realMatrixPart(element: EquationElement) -> EquationElement{
    let result = EquationElement()
    result.type = .MATRIX
    
    for(var index = 0; index < element.matrixComplex.count; index += 2){
        result.matrix.append(element.matrixComplex[index])
    }
    return result
}

func complexMatrixPart(element: EquationElement) -> EquationElement{
    let result = EquationElement()
    result.type = .MATRIX
    
    for(var index = 1; index < element.matrixComplex.count; index += 2){
        result.matrix.append(element.matrixComplex[index])
    }
    return result
}

func convertElementToComplex(element: EquationElement) -> EquationElement{
    let result = EquationElement()
    if(element.type == .VECTOR){
        result.vectorComplex = [Float](count: element.vector.count * 2, repeatedValue: 0.0)
        for(var index = 0; index < element.vector.count; index += 1){
            result.vectorComplex[index * 2] = element.vector[index]
        }
        result.type = .VECTORCOMPLEX
    }
    else if(element.type == .MATRIX){
        result.matrixComplex = [Float](count: element.matrix.count * 2, repeatedValue: 0.0)
        for(var index = 0; index < element.vector.count; index += 1){
            result.matrixComplex[index * 2] = element.matrix[index]
        }
        result.type = .MATRIXCOMPLEX
    }
    else if(element.type == .CONSTANT){
        result.constantComplex = [element.constant, 0.0]
        result.type = .CONSTANTCOMPLEX
    }
    return result
}

func ifElementNotComplexConvertToReal(element: EquationElement) -> EquationElement{//@@@@@@@@@@@@@@@@
    var result = EquationElement()
    if(element.type == .VECTORCOMPLEX){
        var totalComplexValues = 0
        for var index = 1; index < element.vectorComplex.count; index += 2{
            if(element.vectorComplex[index] != 0.0){
                totalComplexValues += 1
            }
        }
        if(totalComplexValues != 0){
            result = element
        }
        else{
            for var index = 0; index < element.vectorComplex.count; index += 2{
                result.vector.append(element.vectorComplex[index])
            }
            result.type = .VECTOR
        }
    }
    else if(element.type == .MATRIXCOMPLEX){
        var totalComplexValues = 0
        for var index = 1; index < element.matrixComplex.count; index += 2{
            if(element.matrixComplex[index] != 0.0){
                totalComplexValues += 1
            }
        }
        if(totalComplexValues != 0){
            result = element
        }
        else{
            for var index = 0; index < element.matrixComplex.count; index += 2{
                result.matrix.append(element.matrixComplex[index])
            }
            result.type = .MATRIX
        }
    }
    else if(element.type == .CONSTANTCOMPLEX){
        if(element.constantComplex[1] != 0.0){
            result = element
        }
        else{
            result.constant = element.constantComplex[0]
            result.type = .CONSTANT
        }
    }
    else{
        result = element
    }
    return result
}

func variableOperatorVariable(firstElement: EquationElement, secondElement: EquationElement, thirdElement: EquationElement) -> EquationElement{
    var result = EquationElement()
    if(secondElement.type == .OPERATION){
        if(secondElement.operation == .ADD){
            if(firstElement.type == .CONSTANT && thirdElement.type == .CONSTANT){
                result.constant = addCC(firstElement.constant, thirdElement.constant)
                result.type = .CONSTANT
            }
            else if(firstElement.type == .CONSTANT && thirdElement.type == .CONSTANTCOMPLEX){
                result.constantComplex = addCCZ(convertElementToComplex(firstElement).constantComplex, thirdElement.constantComplex)
                result.type = .CONSTANTCOMPLEX
            }
            else if(firstElement.type == .CONSTANT && thirdElement.type == .VECTOR){
                result.vector = addVC(thirdElement.vector, firstElement.constant)
                result.type = .VECTOR
            }
            else if(firstElement.type == .CONSTANT && thirdElement.type == .VECTORCOMPLEX){
                result.vectorComplex = addVCZ(thirdElement.vectorComplex, convertElementToComplex(firstElement).constantComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .CONSTANT && thirdElement.type == .MATRIX){
                result.matrix = addMC(thirdElement.matrix, firstElement.constant)
                result.type = .MATRIX
            }
            else if(firstElement.type == .CONSTANT && thirdElement.type == .MATRIXCOMPLEX){
                result.matrixComplex = addMCZ(thirdElement.matrixComplex, convertElementToComplex(firstElement).constantComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .CONSTANTCOMPLEX && thirdElement.type == .CONSTANT){
                result.constantComplex = addCCZ(firstElement.constantComplex, convertElementToComplex(thirdElement).constantComplex)
                result.type = .CONSTANTCOMPLEX
            }
            else if(firstElement.type == .CONSTANTCOMPLEX && thirdElement.type == .CONSTANTCOMPLEX){
                result.constantComplex = addCCZ(firstElement.constantComplex, thirdElement.constantComplex)
                result.type = .CONSTANTCOMPLEX
            }
            else if(firstElement.type == .CONSTANTCOMPLEX && thirdElement.type == .VECTOR){
                result.vectorComplex = addVCZ(convertElementToComplex(thirdElement).vectorComplex, firstElement.constantComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .CONSTANTCOMPLEX && thirdElement.type == .VECTORCOMPLEX){
                result.vectorComplex = addVCZ(thirdElement.vectorComplex, firstElement.constantComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .CONSTANTCOMPLEX && thirdElement.type == .MATRIX){
                result.matrixComplex = addMCZ(convertElementToComplex(thirdElement).matrixComplex, firstElement.constantComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .CONSTANTCOMPLEX && thirdElement.type == .MATRIXCOMPLEX){
                result.matrixComplex = addMCZ(thirdElement.matrixComplex, firstElement.constantComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .VECTOR && thirdElement.type == .CONSTANT){
                result.vector = addVC(firstElement.vector, thirdElement.constant)
                result.type = .VECTOR
            }
            else if(firstElement.type == .VECTOR && thirdElement.type == .CONSTANTCOMPLEX){
                result.vectorComplex = addVCZ(convertElementToComplex(firstElement).vectorComplex, thirdElement.constantComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .VECTOR && thirdElement.type == .VECTOR){
                result.vector = addVV(firstElement.vector, thirdElement.vector)
                result.type = .VECTOR
            }
            else if(firstElement.type == .VECTOR && thirdElement.type == .VECTORCOMPLEX){
                result.vectorComplex = addVVZ(convertElementToComplex(firstElement).vectorComplex, thirdElement.vectorComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .VECTORCOMPLEX && thirdElement.type == .CONSTANT){
                result.vectorComplex = addVCZ(firstElement.vectorComplex, convertElementToComplex(thirdElement).constantComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .VECTORCOMPLEX && thirdElement.type == .CONSTANTCOMPLEX){
                result.vectorComplex = addVCZ(firstElement.vectorComplex, thirdElement.constantComplex)
                result.type = .VECTORCOMPLEX
                result  = ifElementNotComplexConvertToReal(result)
            }
            else if(firstElement.type == .VECTORCOMPLEX && thirdElement.type == .VECTOR){
                result.vectorComplex = addVVZ(firstElement.vectorComplex, convertElementToComplex(thirdElement).vectorComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .VECTORCOMPLEX && thirdElement.type == .VECTORCOMPLEX){
                result.vectorComplex = addVVZ(firstElement.vectorComplex, thirdElement.vectorComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .MATRIX && thirdElement.type == .CONSTANT){
                result.matrix = addMC(firstElement.matrix, thirdElement.constant)
                result.type = .MATRIX
            }
            else if(firstElement.type == .MATRIX && thirdElement.type == .CONSTANTCOMPLEX){
                result.matrixComplex = addMCZ(convertElementToComplex(firstElement).matrixComplex, thirdElement.constantComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .MATRIXCOMPLEX && thirdElement.type == .CONSTANT){
                result.matrixComplex = addMCZ(firstElement.matrixComplex, convertElementToComplex(thirdElement).constantComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .MATRIXCOMPLEX && thirdElement.type == .CONSTANTCOMPLEX){
                result.matrixComplex = addMCZ(firstElement.matrixComplex, thirdElement.constantComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .MATRIX && thirdElement.type == .MATRIX){
                result.matrix = addMM(firstElement.matrix, thirdElement.matrix)
                result.type = .MATRIX
            }
            else if(firstElement.type == .MATRIX && thirdElement.type == .MATRIXCOMPLEX){
                result.matrixComplex = addMMZ(convertElementToComplex(firstElement).matrixComplex, thirdElement.matrixComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .MATRIXCOMPLEX && thirdElement.type == .MATRIX){
                result.matrixComplex = addMMZ(firstElement.matrixComplex, convertElementToComplex(thirdElement).matrixComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .MATRIXCOMPLEX && thirdElement.type == .MATRIXCOMPLEX){
                result.matrixComplex = addMMZ(firstElement.matrixComplex, thirdElement.matrixComplex)
                result.type = .MATRIXCOMPLEX
            }
        }
        else if(secondElement.operation == .SUBTRACT){
            if(firstElement.type == .CONSTANT && thirdElement.type == .CONSTANT){
                result.constant = subtractCC(firstElement.constant, thirdElement.constant)
                result.type = .CONSTANT
            }
            else if(firstElement.type == .CONSTANT && thirdElement.type == .CONSTANTCOMPLEX){
                result.constantComplex = subtractCCZ(convertElementToComplex(firstElement).constantComplex, thirdElement.constantComplex)
                result.type = .CONSTANTCOMPLEX
            }
            else if(firstElement.type == .CONSTANT && thirdElement.type == .VECTOR){
                result.vector = subtractVC(thirdElement.vector, firstElement.constant)
                result.type = .VECTOR
            }
            else if(firstElement.type == .CONSTANT && thirdElement.type == .VECTORCOMPLEX){
                result.vectorComplex = subtractVCZ(thirdElement.vectorComplex, convertElementToComplex(firstElement).constantComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .CONSTANT && thirdElement.type == .MATRIX){
                result.matrix = subtractMC(thirdElement.matrix, firstElement.constant)
                result.type = .MATRIX
            }
            else if(firstElement.type == .CONSTANT && thirdElement.type == .MATRIXCOMPLEX){
                result.matrixComplex = subtractMCZ(thirdElement.matrixComplex, convertElementToComplex(firstElement).constantComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .CONSTANTCOMPLEX && thirdElement.type == .CONSTANT){
                result.constantComplex = subtractCCZ(firstElement.constantComplex, convertElementToComplex(thirdElement).constantComplex)
                result.type = .CONSTANTCOMPLEX
            }
            else if(firstElement.type == .CONSTANTCOMPLEX && thirdElement.type == .CONSTANTCOMPLEX){
                result.constantComplex = subtractCCZ(firstElement.constantComplex, thirdElement.constantComplex)
                result.type = .CONSTANTCOMPLEX
            }
            else if(firstElement.type == .CONSTANTCOMPLEX && thirdElement.type == .VECTOR){
                result.vectorComplex = subtractVCZ(convertElementToComplex(thirdElement).vectorComplex, firstElement.constantComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .CONSTANTCOMPLEX && thirdElement.type == .VECTORCOMPLEX){
                result.vectorComplex = subtractVCZ(thirdElement.vectorComplex, firstElement.constantComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .CONSTANTCOMPLEX && thirdElement.type == .MATRIX){
                result.matrixComplex = subtractMCZ(convertElementToComplex(thirdElement).matrixComplex, firstElement.constantComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .CONSTANTCOMPLEX && thirdElement.type == .MATRIXCOMPLEX){
                result.matrixComplex = subtractMCZ(thirdElement.matrixComplex, firstElement.constantComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .VECTOR && thirdElement.type == .CONSTANT){
                result.vector = subtractVC(firstElement.vector, thirdElement.constant)
                result.type = .VECTOR
            }
            else if(firstElement.type == .VECTOR && thirdElement.type == .CONSTANTCOMPLEX){
                result.vectorComplex = subtractVCZ(convertElementToComplex(firstElement).vectorComplex, thirdElement.constantComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .VECTOR && thirdElement.type == .VECTOR){
                result.vector = subtractVV(firstElement.vector, thirdElement.vector)
                result.type = .VECTOR
            }
            else if(firstElement.type == .VECTOR && thirdElement.type == .VECTORCOMPLEX){
                result.vectorComplex = subtractVVZ(convertElementToComplex(firstElement).vectorComplex, thirdElement.vectorComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .VECTORCOMPLEX && thirdElement.type == .CONSTANT){
                result.vectorComplex = subtractVCZ(firstElement.vectorComplex, convertElementToComplex(thirdElement).constantComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .VECTORCOMPLEX && thirdElement.type == .CONSTANTCOMPLEX){
                result.vectorComplex = subtractVCZ(firstElement.vectorComplex, thirdElement.constantComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .VECTORCOMPLEX && thirdElement.type == .VECTOR){
                result.vectorComplex = subtractVVZ(firstElement.vectorComplex, convertElementToComplex(thirdElement).vectorComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .VECTORCOMPLEX && thirdElement.type == .VECTORCOMPLEX){
                result.vectorComplex = subtractVVZ(firstElement.vectorComplex, thirdElement.vectorComplex)
                result.type = .VECTORCOMPLEX
                result  = ifElementNotComplexConvertToReal(result)
            }
            else if(firstElement.type == .MATRIX && thirdElement.type == .CONSTANT){
                result.matrix = subtractMC(firstElement.matrix, thirdElement.constant)
                result.type = .MATRIX
            }
            else if(firstElement.type == .MATRIX && thirdElement.type == .CONSTANTCOMPLEX){
                result.matrixComplex = subtractMCZ(convertElementToComplex(firstElement).matrixComplex, thirdElement.constantComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .MATRIXCOMPLEX && thirdElement.type == .CONSTANT){
                result.matrixComplex = subtractMCZ(firstElement.matrixComplex, convertElementToComplex(thirdElement).constantComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .MATRIXCOMPLEX && thirdElement.type == .CONSTANTCOMPLEX){
                result.matrixComplex = subtractMCZ(firstElement.matrixComplex, thirdElement.constantComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .MATRIX && thirdElement.type == .MATRIX){
                result.matrix = subtractMM(firstElement.matrix, thirdElement.matrix)
                result.type = .MATRIX
            }
            else if(firstElement.type == .MATRIX && thirdElement.type == .MATRIXCOMPLEX){
                result.matrixComplex = subtractMMZ(convertElementToComplex(firstElement).matrixComplex, thirdElement.matrixComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .MATRIXCOMPLEX && thirdElement.type == .MATRIX){
                result.matrixComplex = subtractMMZ(firstElement.matrixComplex, convertElementToComplex(thirdElement).matrixComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .MATRIXCOMPLEX && thirdElement.type == .MATRIXCOMPLEX){
                result.matrixComplex = subtractMMZ(firstElement.matrixComplex, thirdElement.matrixComplex)
                result.type = .MATRIXCOMPLEX
            }
        }
        else if(secondElement.operation == .MULTIPLY){
            if(firstElement.type == .CONSTANT && thirdElement.type == .CONSTANT){
                result.constant = multiplyCC(firstElement.constant, thirdElement.constant)
                result.type = .CONSTANT
            }
            else if(firstElement.type == .CONSTANT && thirdElement.type == .CONSTANTCOMPLEX){
                result.constantComplex = multiplyCCZ(convertElementToComplex(firstElement).constantComplex, thirdElement.constantComplex)
                result.type = .CONSTANTCOMPLEX
            }
            else if(firstElement.type == .CONSTANT && thirdElement.type == .VECTOR){
                result.vector = multiplyVC(thirdElement.vector, firstElement.constant)
                result.type = .VECTOR
            }
            else if(firstElement.type == .CONSTANT && thirdElement.type == .VECTORCOMPLEX){
                result.vectorComplex = multiplyVCZ(thirdElement.vectorComplex, convertElementToComplex(firstElement).constantComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .CONSTANT && thirdElement.type == .MATRIX){
                result.matrix = multiplyMC(thirdElement.matrix, firstElement.constant)
                result.type = .MATRIX
            }
            else if(firstElement.type == .CONSTANT && thirdElement.type == .MATRIXCOMPLEX){
                result.matrixComplex = multiplyMCZ(thirdElement.matrixComplex, convertElementToComplex(firstElement).constantComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .CONSTANTCOMPLEX && thirdElement.type == .CONSTANT){
                result.constantComplex = multiplyCCZ(firstElement.constantComplex, convertElementToComplex(thirdElement).constantComplex)
                result.type = .CONSTANTCOMPLEX
            }
            else if(firstElement.type == .CONSTANTCOMPLEX && thirdElement.type == .CONSTANTCOMPLEX){
                result.constantComplex = multiplyCCZ(firstElement.constantComplex, thirdElement.constantComplex)
                result.type = .CONSTANTCOMPLEX
            }
            else if(firstElement.type == .CONSTANTCOMPLEX && thirdElement.type == .VECTOR){
                result.vectorComplex = multiplyVCZ(convertElementToComplex(thirdElement).vectorComplex, firstElement.constantComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .CONSTANTCOMPLEX && thirdElement.type == .VECTORCOMPLEX){
                result.vectorComplex = multiplyVCZ(thirdElement.vectorComplex, firstElement.constantComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .CONSTANTCOMPLEX && thirdElement.type == .MATRIX){
                result.matrixComplex = multiplyMCZ(convertElementToComplex(thirdElement).matrixComplex, firstElement.constantComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .CONSTANTCOMPLEX && thirdElement.type == .MATRIXCOMPLEX){
                result.matrixComplex = multiplyMCZ(thirdElement.matrixComplex, firstElement.constantComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .VECTOR && thirdElement.type == .CONSTANT){
                result.vector = multiplyVC(firstElement.vector, thirdElement.constant)
                result.type = .VECTOR
            }
            else if(firstElement.type == .VECTOR && thirdElement.type == .CONSTANTCOMPLEX){
                result.vectorComplex = multiplyVCZ(convertElementToComplex(firstElement).vectorComplex, thirdElement.constantComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .VECTORCOMPLEX && thirdElement.type == .CONSTANT){
                result.vectorComplex = multiplyVCZ(firstElement.vectorComplex, convertElementToComplex(thirdElement).constantComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .VECTORCOMPLEX && thirdElement.type == .CONSTANTCOMPLEX){
                result.vectorComplex = multiplyVCZ(firstElement.vectorComplex, thirdElement.constantComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .MATRIX && thirdElement.type == .CONSTANT){
                result.matrix = multiplyMC(firstElement.matrix, thirdElement.constant)
                result.type = .MATRIX
            }
            else if(firstElement.type == .MATRIX && thirdElement.type == .CONSTANTCOMPLEX){
                result.matrixComplex = multiplyMCZ(convertElementToComplex(firstElement).matrixComplex, thirdElement.constantComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .MATRIX && thirdElement.type == .VECTOR){
                result.vector = multiplyMV(firstElement.matrix, thirdElement.vector)
                result.type = .VECTOR
            }
            else if(firstElement.type == .MATRIX && thirdElement.type == .VECTORCOMPLEX){
                result.vectorComplex = multiplyMVZ(convertElementToComplex(firstElement).matrixComplex, thirdElement.vectorComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .MATRIX && thirdElement.type == .MATRIX){
                result.matrix = multiplyMM(firstElement.matrix, thirdElement.matrix)
                result.type = .MATRIX
            }
            else if(firstElement.type == .MATRIX && thirdElement.type == .MATRIXCOMPLEX){
                result.matrixComplex = multiplyMMZ(convertElementToComplex(firstElement).matrixComplex, thirdElement.matrixComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .MATRIXCOMPLEX && thirdElement.type == .CONSTANT){
                result.matrixComplex = multiplyMCZ(firstElement.matrixComplex, convertElementToComplex(thirdElement).constantComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .MATRIXCOMPLEX && thirdElement.type == .CONSTANTCOMPLEX){
                result.matrixComplex = multiplyMCZ(firstElement.matrixComplex, thirdElement.constantComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .MATRIXCOMPLEX && thirdElement.type == .VECTOR){
                result.vectorComplex = multiplyMVZ(firstElement.matrixComplex, convertElementToComplex(thirdElement).vectorComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .MATRIXCOMPLEX && thirdElement.type == .VECTORCOMPLEX){
                result.vectorComplex = multiplyMVZ(firstElement.matrixComplex, thirdElement.vectorComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .MATRIXCOMPLEX && thirdElement.type == .MATRIX){
                result.matrixComplex = multiplyMMZ(firstElement.matrixComplex, convertElementToComplex(thirdElement).matrixComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .MATRIXCOMPLEX && thirdElement.type == .MATRIXCOMPLEX){
                result.matrixComplex = multiplyMMZ(firstElement.matrixComplex, thirdElement.matrixComplex)
                result.type = .MATRIXCOMPLEX
            }
        }
        else if(secondElement.operation == .DIVIDE){
            if(firstElement.type == .CONSTANT && thirdElement.type == .CONSTANT){
                result.constant = divideCC(firstElement.constant, thirdElement.constant)
                result.type = .CONSTANT
            }
            else if(firstElement.type == .CONSTANT && thirdElement.type == .CONSTANTCOMPLEX){
                result.constantComplex = divideCCZ(convertElementToComplex(firstElement).constantComplex, thirdElement.constantComplex)
                result.type = .CONSTANTCOMPLEX
            }
            else if(firstElement.type == .CONSTANT && thirdElement.type == .VECTOR){
                result.vector = divideVC(thirdElement.vector, firstElement.constant)
                result.type = .VECTOR
            }
            else if(firstElement.type == .CONSTANT && thirdElement.type == .VECTORCOMPLEX){
                result.vectorComplex = divideVCZ(thirdElement.vectorComplex, convertElementToComplex(firstElement).constantComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .CONSTANT && thirdElement.type == .MATRIX){
                result.matrix = divideMC(thirdElement.matrix, firstElement.constant)
                result.type = .MATRIX
            }
            else if(firstElement.type == .CONSTANT && thirdElement.type == .MATRIXCOMPLEX){
                result.matrixComplex = divideMCZ(thirdElement.matrixComplex, convertElementToComplex(firstElement).constantComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .CONSTANTCOMPLEX && thirdElement.type == .CONSTANT){
                result.constantComplex = divideCCZ(firstElement.constantComplex, convertElementToComplex(thirdElement).constantComplex)
                result.type = .CONSTANTCOMPLEX
            }
            else if(firstElement.type == .CONSTANTCOMPLEX && thirdElement.type == .CONSTANTCOMPLEX){
                result.constantComplex = divideCCZ(firstElement.constantComplex, thirdElement.constantComplex)
                result.type = .CONSTANTCOMPLEX
            }
            else if(firstElement.type == .CONSTANTCOMPLEX && thirdElement.type == .VECTOR){
                result.vectorComplex = divideVCZ(convertElementToComplex(thirdElement).vectorComplex, firstElement.constantComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .CONSTANTCOMPLEX && thirdElement.type == .VECTORCOMPLEX){
                result.vectorComplex = divideVCZ(thirdElement.vectorComplex, firstElement.constantComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .CONSTANTCOMPLEX && thirdElement.type == .MATRIX){
                result.matrixComplex = divideMCZ(convertElementToComplex(thirdElement).matrixComplex, firstElement.constantComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .CONSTANTCOMPLEX && thirdElement.type == .MATRIXCOMPLEX){
                result.matrixComplex = divideMCZ(thirdElement.matrixComplex, firstElement.constantComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .VECTOR && thirdElement.type == .CONSTANT){
                result.vector = divideVC(firstElement.vector, thirdElement.constant)
                result.type = .VECTOR
            }
            else if(firstElement.type == .VECTOR && thirdElement.type == .CONSTANTCOMPLEX){
                result.vectorComplex = divideVCZ(convertElementToComplex(firstElement).vectorComplex, thirdElement.constantComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .VECTORCOMPLEX && thirdElement.type == .CONSTANT){
                result.vectorComplex = divideVCZ(firstElement.vectorComplex, convertElementToComplex(thirdElement).constantComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .VECTORCOMPLEX && thirdElement.type == .CONSTANTCOMPLEX){
                result.vectorComplex = divideVCZ(firstElement.vectorComplex, thirdElement.constantComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .MATRIX && thirdElement.type == .CONSTANT){
                result.matrix = divideMC(firstElement.matrix, thirdElement.constant)
                result.type = .MATRIX
            }
            else if(firstElement.type == .MATRIX && thirdElement.type == .CONSTANTCOMPLEX){
                result.matrixComplex = divideMCZ(convertElementToComplex(firstElement).matrixComplex, thirdElement.constantComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .MATRIXCOMPLEX && thirdElement.type == .CONSTANT){
                result.matrixComplex = divideMCZ(firstElement.matrixComplex, convertElementToComplex(thirdElement).constantComplex)
                result.type = .MATRIXCOMPLEX
            }
            else if(firstElement.type == .MATRIXCOMPLEX && thirdElement.type == .CONSTANTCOMPLEX){
                result.matrixComplex = divideMCZ(firstElement.matrixComplex, thirdElement.constantComplex)
                result.type = .MATRIXCOMPLEX
            }
        }
    }
    result  = ifElementNotComplexConvertToReal(result)
    return result
}

func operatorVariable(firstElement: EquationElement, secondElement: EquationElement) -> EquationElement{
    var result = EquationElement()
    if(firstElement.type == .OPERATION){
        if(firstElement.operation == .NEGATE){
            if(secondElement.type == .CONSTANT){
                result.constant = negateC(secondElement.constant)
                result.type = .CONSTANT
            }
            else if(secondElement.type == .CONSTANTCOMPLEX){
                result.constantComplex = negateCZ(secondElement.constantComplex)
                result.type = .CONSTANTCOMPLEX
            }
            else if(secondElement.type == .VECTOR){
                result.vector = negateV(secondElement.vector)
                result.type = .VECTOR
            }
            else if(secondElement.type == .VECTORCOMPLEX){
                result.vectorComplex = negateVZ(secondElement.vectorComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(secondElement.type == .MATRIX){
                result.matrix = negateM(secondElement.matrix)
                result.type = .MATRIX
            }
            else if(secondElement.type == .MATRIXCOMPLEX){
                result.matrixComplex = negateMZ(secondElement.matrixComplex)
                result.type = .MATRIXCOMPLEX
            }
        }
        else if(firstElement.operation == .SQUAREROOT){
            if(secondElement.type == .CONSTANT){
                result.constant = squareRootC(secondElement.constant)
                result.type = .CONSTANT
            }
        }
    }
    result  = ifElementNotComplexConvertToReal(result)
    return result
}

func variableOperator(firstElement: EquationElement, secondElement: EquationElement) -> EquationElement{
    var result = EquationElement()
    if(secondElement.type == .OPERATION){
        if(secondElement.operation == .TRANSPOSE){
            if(firstElement.type == .MATRIX){
                result.matrix = transposeM(firstElement.matrix)
                result.type = .MATRIX
            }
            else if(firstElement.type == .MATRIXCOMPLEX){
                result.matrixComplex = transposeMZ(firstElement.matrixComplex)
                result.type = .MATRIXCOMPLEX
            }
        }
        else if(secondElement.operation == .COMPLEXCONJUGATE){
            if(firstElement.type == .CONSTANTCOMPLEX){
                result.constantComplex = complexConjugateCZ(firstElement.constantComplex)
                result.type = .CONSTANTCOMPLEX
            }
            else if(firstElement.type == .VECTORCOMPLEX){
                result.vectorComplex = complexConjugateVZ(firstElement.vectorComplex)
                result.type = .VECTORCOMPLEX
            }
            else if(firstElement.type == .MATRIXCOMPLEX){
                result.matrixComplex = complexConjugateMZ(firstElement.matrixComplex)
                result.type = .MATRIXCOMPLEX
            }
        }
        else if(secondElement.operation == .DAGGER){
            if(firstElement.type == .MATRIX){
                result.matrix = daggerM(firstElement.matrix)
                result.type = .MATRIX
            }
            else if(firstElement.type == .MATRIXCOMPLEX){
                result.matrixComplex = daggerMZ(firstElement.matrixComplex)
                result.type = .MATRIXCOMPLEX
            }
        }
        else if(secondElement.operation == .SQUARED){
            if(firstElement.type == .CONSTANT){
                result.constant = squaredC(firstElement.constant)
                result.type = .CONSTANT
                result  = ifElementNotComplexConvertToReal(result)
            }
            else if(firstElement.type == .MATRIX){
                result.matrix = squaredM(firstElement.matrix)
                result.type = .MATRIX
            }
            else if(firstElement.type == .CONSTANTCOMPLEX){
                result.constantComplex = squaredCZ(firstElement.constantComplex)
                result.type = .CONSTANTCOMPLEX
            }
            else if(firstElement.type == .MATRIXCOMPLEX){
                result.matrixComplex = squaredMZ(firstElement.matrixComplex)
                result.type = .MATRIXCOMPLEX
            }
        }
        else if(secondElement.operation == .CUBED){
            if(firstElement.type == .CONSTANT){
                result.constant = cubedC(firstElement.constant)
                result.type = .CONSTANT
            }
            else if(firstElement.type == .MATRIX){
                result.matrix = cubedM(firstElement.matrix)
                result.type = .MATRIX
            }
            else if(firstElement.type == .CONSTANTCOMPLEX){
                result.constantComplex = cubedCZ(firstElement.constantComplex)
                result.type = .CONSTANTCOMPLEX
            }
            else if(firstElement.type == .MATRIXCOMPLEX){
                result.matrixComplex = cubedMZ(firstElement.matrixComplex)
                result.type = .MATRIXCOMPLEX
            }
        }
    }
    result  = ifElementNotComplexConvertToReal(result)
    return result
}

func implimentSquareRoot(elementArray: [EquationElement]) -> [EquationElement]{
    var totalSquareRoot = 0
    for element in elementArray{
        if(element.operation == .SQUAREROOT){
            totalSquareRoot += 1
        }
    }
    
    var resultElementArray: [EquationElement] = elementArray
    
    while(totalSquareRoot > 0){
        var newElementArray = [EquationElement](count: resultElementArray.count - 1, repeatedValue: EquationElement())
        var newElementArrayIndex = 0
        
        var squareRootIndex: Int = -1
        for var index = 0; index < resultElementArray.count; index += 1{
            if(resultElementArray[index].operation == .SQUAREROOT){
                squareRootIndex = index
            }
        }
        
        for var index = 0; index < squareRootIndex; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        newElementArray[newElementArrayIndex] = operatorVariable(resultElementArray[squareRootIndex], resultElementArray[squareRootIndex+1])
        newElementArrayIndex += 1
        
        for var index = squareRootIndex + 2; index < resultElementArray.count; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        resultElementArray = newElementArray
        totalSquareRoot -= 1
    }
    return resultElementArray
}

func implimentNegate(elementArray: [EquationElement]) -> [EquationElement]{
    var totalNegate = 0
    for element in elementArray{
        if(element.operation == .NEGATE){
            totalNegate += 1
        }
    }
    
    var resultElementArray: [EquationElement] = elementArray
    
    while(totalNegate > 0){
        var newElementArray = [EquationElement](count: resultElementArray.count - 1, repeatedValue: EquationElement())
        var newElementArrayIndex = 0

        var negateIndex: Int = -1
        for var index = 0; index < resultElementArray.count; index += 1{
            if(resultElementArray[index].operation == .NEGATE){
                negateIndex = index
            }
        }
        
        for var index = 0; index < negateIndex; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        newElementArray[newElementArrayIndex] = operatorVariable(resultElementArray[negateIndex], resultElementArray[negateIndex+1])
        newElementArrayIndex += 1
        
        for var index = negateIndex + 2; index < resultElementArray.count; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        resultElementArray = newElementArray
        totalNegate -= 1
    }
    return resultElementArray
}

func implimentTranspose(elementArray: [EquationElement]) -> [EquationElement]{
    var totalTranspose = 0
    for element in elementArray{
        if(element.operation == .TRANSPOSE){
            totalTranspose += 1
        }
    }
    
    var resultElementArray: [EquationElement] = elementArray
    
    while(totalTranspose > 0){
        var newElementArray = [EquationElement](count: resultElementArray.count - 1, repeatedValue: EquationElement())
        var newElementArrayIndex = 0

        var transposeIndex: Int = -1
        for var index = 0; index < resultElementArray.count; index += 1{
            if(resultElementArray[index].operation == .TRANSPOSE){
                transposeIndex = index
            }
        }
        
        for var index = 0; index < transposeIndex - 1; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        newElementArray[newElementArrayIndex] = variableOperator(resultElementArray[transposeIndex-1], resultElementArray[transposeIndex])
        newElementArrayIndex += 1
        
        for var index = transposeIndex + 1; index < resultElementArray.count; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        resultElementArray = newElementArray
        totalTranspose -= 1
    }
    return resultElementArray
}

func implimentComplexConjugate(elementArray: [EquationElement]) -> [EquationElement]{
    var totalComplexConjugate = 0
    for element in elementArray{
        if(element.operation == .COMPLEXCONJUGATE){
            totalComplexConjugate += 1
        }
    }
    
    var resultElementArray: [EquationElement] = elementArray
    
    while(totalComplexConjugate > 0){
        var newElementArray = [EquationElement](count: resultElementArray.count - 1, repeatedValue: EquationElement())
        var newElementArrayIndex = 0
        
        var complexConjugateIndex: Int = -1
        for var index = 0; index < resultElementArray.count; index += 1{
            if(resultElementArray[index].operation == .COMPLEXCONJUGATE){
                complexConjugateIndex = index
            }
        }
        
        for var index = 0; index < complexConjugateIndex - 1; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        newElementArray[newElementArrayIndex] = variableOperator(resultElementArray[complexConjugateIndex-1], resultElementArray[complexConjugateIndex])
        newElementArrayIndex += 1
        
        for var index = complexConjugateIndex + 1; index < resultElementArray.count; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        resultElementArray = newElementArray
        totalComplexConjugate -= 1
    }
    return resultElementArray
}

func implimentDagger(elementArray: [EquationElement]) -> [EquationElement]{
    var totalDagger = 0
    for element in elementArray{
        if(element.operation == .DAGGER){
            totalDagger += 1
        }
    }
    
    var resultElementArray: [EquationElement] = elementArray
    
    while(totalDagger > 0){
        var newElementArray = [EquationElement](count: resultElementArray.count - 1, repeatedValue: EquationElement())
        var newElementArrayIndex = 0

        var daggerIndex: Int = -1
        for var index = 0; index < resultElementArray.count; index += 1{
            if(resultElementArray[index].operation == .DAGGER){
                daggerIndex = index
            }
        }
        
        for var index = 0; index < daggerIndex - 1; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        newElementArray[newElementArrayIndex] = variableOperator(resultElementArray[daggerIndex-1], resultElementArray[daggerIndex])
        newElementArrayIndex += 1
        
        for var index = daggerIndex + 1; index < resultElementArray.count; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        resultElementArray = newElementArray
        totalDagger -= 1
    }
    return resultElementArray
}

func implimentSquared(elementArray: [EquationElement]) -> [EquationElement]{
    var totalSquared = 0
    for element in elementArray{
        if(element.operation == .SQUARED){
            totalSquared += 1
        }
    }
    
    var resultElementArray: [EquationElement] = elementArray
    
    while(totalSquared > 0){
        var newElementArray = [EquationElement](count: resultElementArray.count - 1, repeatedValue: EquationElement())
        var newElementArrayIndex = 0
        
        var squaredIndex: Int = -1
        for var index = 0; index < resultElementArray.count; index += 1{
            if(resultElementArray[index].operation == .SQUARED){
                squaredIndex = index
            }
        }
        
        for var index = 0; index < squaredIndex - 1; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        newElementArray[newElementArrayIndex] = variableOperator(resultElementArray[squaredIndex-1], resultElementArray[squaredIndex])
        newElementArrayIndex += 1
        
        for var index = squaredIndex + 1; index < resultElementArray.count; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        resultElementArray = newElementArray
        totalSquared -= 1
    }
    return resultElementArray
}

func implimentCubed(elementArray: [EquationElement]) -> [EquationElement]{
    var totalCubed = 0
    for element in elementArray{
        if(element.operation == .CUBED){
            totalCubed += 1
        }
    }
    
    var resultElementArray: [EquationElement] = elementArray
    
    while(totalCubed > 0){
        var newElementArray = [EquationElement](count: resultElementArray.count - 1, repeatedValue: EquationElement())
        var newElementArrayIndex = 0
        
        var cubedIndex: Int = -1
        for var index = 0; index < resultElementArray.count; index += 1{
            if(resultElementArray[index].operation == .CUBED){
                cubedIndex = index
            }
        }
        
        for var index = 0; index < cubedIndex - 1; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        newElementArray[newElementArrayIndex] = variableOperator(resultElementArray[cubedIndex-1], resultElementArray[cubedIndex])
        newElementArrayIndex += 1
        
        for var index = cubedIndex + 1; index < resultElementArray.count; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        resultElementArray = newElementArray
        totalCubed -= 1
    }
    return resultElementArray
}

func power(element: EquationElement, powerValue: EquationElement) -> EquationElement{
    var finalValue: EquationElement = element
    if(powerValue.constant == 0 && element.type == .CONSTANT){
        finalValue = EquationElement(1.0)
    }
    else if(powerValue.constant == 0 && element.type == .CONSTANTCOMPLEX){
        finalValue = EquationElement(1.0)
    }
    else if(powerValue.constant == 0 && element.type == .MATRIX){
        finalValue = identityMatrix(Int(sqrt(Float(element.matrix.count))))
    }
    else if(powerValue.constant == 0 && element.type == .MATRIXCOMPLEX){
        finalValue = identityMatrix(Int(sqrt(Float(element.matrixComplex.count/2))))
    }
    else if(powerValue.constant > 1){
        while(powerValue.constant > 1){
            finalValue = implimentElementArraySimplification([finalValue, EquationElement(.MULTIPLY), element])
            powerValue.constant -= 1
        }
    }
    return finalValue
}

func implimentPower(elementArray: [EquationElement]) -> [EquationElement]{
    var totalPower = 0
    for element in elementArray{
        if(element.operation == .POWER){
            totalPower += 1
        }
    }
    
    var resultElementArray: [EquationElement] = elementArray
    
    while(totalPower > 0){
        var newElementArray = [EquationElement](count: resultElementArray.count - 2, repeatedValue: EquationElement())
        var newElementArrayIndex = 0

        var powerIndex: Int = -1
        for var index = 0; index < resultElementArray.count; index += 1{
            if(resultElementArray[index].operation == .POWER){
                powerIndex = index
            }
        }
        
        for var index = 0; index < powerIndex - 1; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        newElementArray[newElementArrayIndex] = power(resultElementArray[powerIndex-1], resultElementArray[powerIndex+1])
        newElementArrayIndex += 1
        
        for var index = powerIndex + 2; index < resultElementArray.count; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        resultElementArray = newElementArray
        totalPower -= 1
    }
    return resultElementArray
}

func factorial(element: EquationElement) -> EquationElement{
    var total: Float = 1.0
    if(element.constant != 0){
        for var counter = element.constant; counter > 1; counter -= 1{
            total *= counter
        }
    }
    return EquationElement(total)
}

func implimentFactorial(elementArray: [EquationElement]) -> [EquationElement]{
    var totalFactorial = 0
    for element in elementArray{
        if(element.operation == .FACTORIAL){
            totalFactorial += 1
        }
    }
    
    var resultElementArray: [EquationElement] = elementArray
    
    while(totalFactorial > 0){
        var newElementArray = [EquationElement](count: resultElementArray.count - 1, repeatedValue: EquationElement())
        var newElementArrayIndex = 0
        
        var factorialIndex: Int = -1
        for var index = 0; index < resultElementArray.count; index += 1{
            if(resultElementArray[index].operation == .FACTORIAL){
                factorialIndex = index
            }
        }
        
        for var index = 0; index < factorialIndex - 1; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        newElementArray[newElementArrayIndex] = factorial(resultElementArray[factorialIndex-1])
        newElementArrayIndex += 1
        
        for var index = factorialIndex + 1; index < resultElementArray.count; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        resultElementArray = newElementArray
        totalFactorial -= 1
    }
    return resultElementArray
}

func implimentCos(elementArray: [EquationElement]) -> [EquationElement]{
    var totalCos = 0
    for element in elementArray{
        if(element.operation == .COS){
            totalCos += 1
        }
    }
    
    var resultElementArray: [EquationElement] = elementArray
    
    while(totalCos > 0){
        var newElementArray = [EquationElement](count: resultElementArray.count - 1, repeatedValue: EquationElement())
        var newElementArrayIndex = 0
        
        var cosIndex: Int = -1
        for var index = 0; index < resultElementArray.count; index += 1{
            if(resultElementArray[index].operation == .COS){
                cosIndex = index
            }
        }
        
        for var index = 0; index < cosIndex; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        newElementArray[newElementArrayIndex] = EquationElement(cos(resultElementArray[cosIndex+1].constant))
        newElementArrayIndex += 1
        
        for var index = cosIndex + 2; index < resultElementArray.count; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        resultElementArray = newElementArray
        totalCos -= 1
    }
    return resultElementArray
}

func implimentSin(elementArray: [EquationElement]) -> [EquationElement]{
    var totalSin = 0
    for element in elementArray{
        if(element.operation == .SIN){
            totalSin += 1
        }
    }
    
    var resultElementArray: [EquationElement] = elementArray
    
    while(totalSin > 0){
        var newElementArray = [EquationElement](count: resultElementArray.count - 1, repeatedValue: EquationElement())
        var newElementArrayIndex = 0

        var sinIndex: Int = -1
        for var index = 0; index < resultElementArray.count; index += 1{
            if(resultElementArray[index].operation == .SIN){
                sinIndex = index
            }
        }
        
        for var index = 0; index < sinIndex; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        newElementArray[newElementArrayIndex] = EquationElement(sin(resultElementArray[sinIndex+1].constant))
        newElementArrayIndex += 1
        
        for var index = sinIndex + 2; index < resultElementArray.count; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        resultElementArray = newElementArray
        totalSin -= 1
    }
    return resultElementArray
}

func implimentExp(elementArray: [EquationElement]) -> [EquationElement]{
    var totalExp = 0
    for element in elementArray{
        if(element.operation == .EXP){
            totalExp += 1
        }
    }
    
    var resultElementArray: [EquationElement] = elementArray
    
    while(totalExp > 0){
        var newElementArray = [EquationElement](count: resultElementArray.count - 1, repeatedValue: EquationElement())
        var newElementArrayIndex = 0

        var expIndex: Int = -1
        for var index = 0; index < resultElementArray.count; index += 1{
            if(resultElementArray[index].operation == .EXP){
                expIndex = index
            }
        }
        
        for var index = 0; index < expIndex; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        newElementArray[newElementArrayIndex] = EquationElement(exp(resultElementArray[expIndex+1].constant))
        newElementArrayIndex += 1
        
        for var index = expIndex + 2; index < resultElementArray.count; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        resultElementArray = newElementArray
        totalExp -= 1
    }
    return resultElementArray
}

func implimentCosq(elementArray: [EquationElement]) -> [EquationElement]{
    var totalCosq = 0
    for element in elementArray{
        if(element.operation == .COSQ){
            totalCosq += 1
        }
    }
    
    var resultElementArray: [EquationElement] = elementArray
    
    while(totalCosq > 0){
        var newElementArray = [EquationElement](count: resultElementArray.count - 3, repeatedValue: EquationElement())
        var newElementArrayIndex = 0

        var cosqIndex: Int = -1
        for var index = 0; index < resultElementArray.count; index += 1{
            if(resultElementArray[index].operation == .COSQ){
                cosqIndex = index
            }
        }
        
        for var index = 0; index < cosqIndex; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        newElementArray[newElementArrayIndex] = cosq(Int(resultElementArray[cosqIndex+1].constant), resultElementArray[cosqIndex+2].constant, resultElementArray[cosqIndex+3].constant)
        newElementArrayIndex += 1
        
        for var index = cosqIndex + 4; index < resultElementArray.count; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        resultElementArray = newElementArray
        totalCosq -= 1
    }
    return resultElementArray
}


func implimentDivide(elementArray: [EquationElement]) -> [EquationElement]{
    var totalDivide = 0
    for element in elementArray{
        if(element.operation == .DIVIDE){
            totalDivide += 1
        }
    }
    
    var resultElementArray: [EquationElement] = elementArray
    
    while(totalDivide > 0){
        var newElementArray = [EquationElement](count: resultElementArray.count - 2, repeatedValue: EquationElement())
        var newElementArrayIndex = 0

        var divideIndex: Int = -1
        for var index = 0; index < resultElementArray.count; index += 1{
            if(resultElementArray[index].operation == .DIVIDE){
                divideIndex = index
            }
        }
        
        for var index = 0; index < divideIndex - 1; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        newElementArray[newElementArrayIndex] = variableOperatorVariable(resultElementArray[divideIndex-1], resultElementArray[divideIndex], resultElementArray[divideIndex+1])
        newElementArrayIndex += 1
        
        for var index = divideIndex + 2; index < resultElementArray.count; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        resultElementArray = newElementArray
        totalDivide -= 1
    }
    return resultElementArray
}

func implimentMultiply(elementArray: [EquationElement]) -> [EquationElement]{
    var totalMultiply = 0
    for element in elementArray{
        if(element.operation == .MULTIPLY){
            totalMultiply += 1
        }
    }
    
    var resultElementArray: [EquationElement] = elementArray
    
    while(totalMultiply > 0){
        var newElementArray = [EquationElement](count: resultElementArray.count - 2, repeatedValue: EquationElement())
        var newElementArrayIndex = 0
        
        var multiplyIndex: Int = -1
        for var index = 0; index < resultElementArray.count; index += 1{
            if(resultElementArray[index].operation == .MULTIPLY){
                multiplyIndex = index
            }
        }
        
        for var index = 0; index < multiplyIndex - 1; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        newElementArray[newElementArrayIndex] = variableOperatorVariable(resultElementArray[multiplyIndex-1], resultElementArray[multiplyIndex], resultElementArray[multiplyIndex+1])
        newElementArrayIndex += 1
        
        for var index = multiplyIndex + 2; index < resultElementArray.count; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        resultElementArray = newElementArray
        totalMultiply -= 1
    }
    return resultElementArray
}

func implimentAdd(elementArray: [EquationElement]) -> [EquationElement]{
    var totalAdd = 0
    for element in elementArray{
        if(element.operation == .ADD){
            totalAdd += 1
        }
    }
    
    var resultElementArray: [EquationElement] = elementArray
    
    while(totalAdd > 0){
        var newElementArray = [EquationElement](count: resultElementArray.count - 2, repeatedValue: EquationElement())
        var newElementArrayIndex = 0
        
        var addIndex: Int = -1
        for var index = 0; index < resultElementArray.count; index += 1{
            if(resultElementArray[index].operation == .ADD){
                addIndex = index
            }
        }
        
        for var index = 0; index < addIndex - 1; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        newElementArray[newElementArrayIndex] = variableOperatorVariable(resultElementArray[addIndex-1], resultElementArray[addIndex], resultElementArray[addIndex+1])
        newElementArrayIndex += 1
        
        for var index = addIndex + 2; index < resultElementArray.count; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        resultElementArray = newElementArray
        totalAdd -= 1
    }
    return resultElementArray
}

func implimentSubtract(elementArray: [EquationElement]) -> [EquationElement]{
    var totalSubtract = 0
    for element in elementArray{
        if(element.operation == .SUBTRACT){
            totalSubtract += 1
        }
    }
    
    var resultElementArray: [EquationElement] = elementArray
    
    while(totalSubtract > 0){
        var newElementArray = [EquationElement](count: resultElementArray.count - 2, repeatedValue: EquationElement())
        var newElementArrayIndex = 0

        var subtractIndex: Int = -1
        for var index = 0; index < resultElementArray.count; index += 1{
            if(resultElementArray[index].operation == .SUBTRACT){
                subtractIndex = index
            }
        }
        
        for var index = 0; index < subtractIndex - 1; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        newElementArray[newElementArrayIndex] = variableOperatorVariable(resultElementArray[subtractIndex-1], resultElementArray[subtractIndex], resultElementArray[subtractIndex+1])
        newElementArrayIndex += 1
        
        for var index = subtractIndex + 2; index < resultElementArray.count; index += 1{
            newElementArray[newElementArrayIndex] = resultElementArray[index]
            newElementArrayIndex += 1
        }
        
        resultElementArray = newElementArray
        totalSubtract -= 1
    }
    return resultElementArray
}

func implimentElementArrayOperationSolve(elementArray: [EquationElement]) -> EquationElement{
    let total = implimentSubtract(implimentAdd(implimentMultiply(implimentDivide(implimentDagger(implimentComplexConjugate(implimentTranspose(implimentNegate(implimentCubed(implimentSquared(implimentPower(implimentFactorial(implimentSquareRoot(implimentCosq(implimentCos(implimentSin(implimentExp(elementArray)))))))))))))))))
    return total[0]
}

func implimentElementArrayBracketSolve(elementArray: [EquationElement]) -> [EquationElement]{
    var newElementArray: [EquationElement] = []
    var leftBracketIndex: Int = -1
    var rightBracketIndex: Int = -1
    for var index = 0; index < elementArray.count; index += 1{
        if(elementArray[index].braket == .LEFT){
            leftBracketIndex = index
        }
        else if(elementArray[index].braket == .RIGHT){
            rightBracketIndex = index
            break
        }
    }
    
    for var index = 0; index < leftBracketIndex; index += 1{
        newElementArray.append(elementArray[index])
    }
    
    var tempArray: [EquationElement] = []
    for var index = leftBracketIndex + 1; index < rightBracketIndex; index += 1{
        tempArray.append(elementArray[index])
    }
    newElementArray.append(implimentElementArrayOperationSolve(tempArray))
    
    for var index = rightBracketIndex + 1; index < elementArray.count; index += 1{
        newElementArray.append(elementArray[index])
    }
    return newElementArray
}

func implimentElementArraySimplification(elementArray: [EquationElement]) -> EquationElement{
    var finalElementArray: [EquationElement] = elementArray
    var numberOfLeftBrackets = 0
    for element in elementArray{
        if(element.braket == .LEFT){
            numberOfLeftBrackets += 1
        }
    }
    while(numberOfLeftBrackets != 0){
        //print(numberOfLeftBrackets)
        finalElementArray = implimentElementArrayBracketSolve(finalElementArray)
        numberOfLeftBrackets -= 1
    }
    
    return implimentElementArrayOperationSolve(finalElementArray)
}

func reduceMatrixSizeBy1(element: EquationElement) -> EquationElement{
    var counter = 0.0
    for _ in element.matrix{
        counter += 1
    }
    let size = Int(sqrt(counter))
    let newLength = (size - 1)*(size - 1)
    var currentLength = 0
    let newMatrix = EquationElement()
    newMatrix.type = .MATRIX
    for(var value = 0; value < element.matrix.count; value += 1){
        if((value + 1) % size == 0){
            continue
        }
        else{
            if (currentLength < newLength){
                newMatrix.matrix.append(element.matrix[value])
                currentLength += 1
            }
        }
    }
    return newMatrix
}

/*
func eigenValuesVectors(element: EquationElement){
    var matrix = element.matrix
    
    var N = __CLPK_integer(sqrt(Double(matrix.count)))
    var workspace = [Float](count: Int(N), repeatedValue: 0.0)
    var error : __CLPK_integer = 0
    var lwork = __CLPK_integer(-1)
    // Real parts of eigenvalues
    var wr = [Float](count: Int(N), repeatedValue: 0)
    // Imaginary parts of eigenvalues
    var wi = [Float](count: Int(N), repeatedValue: 0)
    // Left eigenvectors
    var vl = [__CLPK_real](count: Int(N*N), repeatedValue: 0)
    // Right eigenvectors
    var vr = [__CLPK_real](count: Int(N*N), repeatedValue: 0)
    
    var workspaceQuery: Float = 0.0
    sgeev_(UnsafeMutablePointer(("V" as NSString).UTF8String), UnsafeMutablePointer(("V" as NSString).UTF8String), &N, &matrix, &N, &wr, &wi, &vl, &N, &vr, &N, &workspaceQuery, &lwork, &error)
    
    // size workspace per the results of the query:
    workspace = [Float](count: Int(workspaceQuery), repeatedValue: 0.0)
    lwork = __CLPK_integer(workspaceQuery)
    
    sgeev_(UnsafeMutablePointer(("V" as NSString).UTF8String), UnsafeMutablePointer(("V" as NSString).UTF8String), &N, &matrix, &N, &wr, &wi, &vl, &N, &vr, &N, &workspace, &lwork, &error)
    
    //print("\(wr), \(vl), \(vr)")
    //print(wr)
}*/

func eigenValues(element: EquationElement) -> [Float]{
    var matrix = element.matrix
    
    var N = __CLPK_integer(sqrt(Double(matrix.count)))
    var workspace = [Float](count: Int(N), repeatedValue: 0.0)
    var error : __CLPK_integer = 0
    var lwork = __CLPK_integer(-1)
    // Real parts of eigenvalues
    var wr = [Float](count: Int(N), repeatedValue: 0)
    // Imaginary parts of eigenvalues
    var wi = [Float](count: Int(N), repeatedValue: 0)
    // Left eigenvectors
    var vl = [__CLPK_real](count: Int(N*N), repeatedValue: 0)
    // Right eigenvectors
    var vr = [__CLPK_real](count: Int(N*N), repeatedValue: 0)
    
    var workspaceQuery: Float = 0.0
    sgeev_(UnsafeMutablePointer(("V" as NSString).UTF8String), UnsafeMutablePointer(("V" as NSString).UTF8String), &N, &matrix, &N, &wr, &wi, &vl, &N, &vr, &N, &workspaceQuery, &lwork, &error)
    
    // size workspace per the results of the query:
    workspace = [Float](count: Int(workspaceQuery), repeatedValue: 0.0)
    lwork = __CLPK_integer(workspaceQuery)
    
    sgeev_(UnsafeMutablePointer(("V" as NSString).UTF8String), UnsafeMutablePointer(("V" as NSString).UTF8String), &N, &matrix, &N, &wr, &wi, &vl, &N, &vr, &N, &workspace, &lwork, &error)
    
    return wr
}



let mass = EquationElement(6.62607e-34)
//let mass = EquationElement(1.0)
let angularMomentum = EquationElement(1.0)
let h = EquationElement(6.62607e-34)
let pi = EquationElement(3.14159265359)
let hBar = implimentElementArraySimplification([h, EquationElement(.DIVIDE), pi])
//let hBar = EquationElement(1.0)
let i = EquationElement(.CONSTANTCOMPLEX, [0.0,1.0])


func annihilationOperator(operatorSize: Int) -> EquationElement{
    let result = EquationElement()
    var tempMatrix = [Float](count: operatorSize * operatorSize, repeatedValue: 0.0)
    var positionIndex = 1
    for(var root = 1; root < operatorSize; root += 1){
        tempMatrix[positionIndex] = Float(sqrt(Double(root)))
        positionIndex += operatorSize + 1
    }
    result.matrix = tempMatrix
    result.type = .MATRIX
    return result
}

func creationOperator(operatorSize: Int) -> EquationElement{
    return implimentElementArraySimplification([annihilationOperator(operatorSize), EquationElement(.DAGGER)])
}

func numberOperator(var operatorSize: Int) -> EquationElement{
    operatorSize += 1
    return reduceMatrixSizeBy1(implimentElementArraySimplification([creationOperator(operatorSize), EquationElement(.MULTIPLY), annihilationOperator(operatorSize)]))
}

func positionOperator(operatorSize: Int) -> EquationElement{
    return implimentElementArraySimplification([EquationElement(.SQUAREROOT), EquationElement(.LEFT), hBar, EquationElement(.DIVIDE), EquationElement(.LEFT), EquationElement(2.0), EquationElement(.MULTIPLY), mass, EquationElement(.MULTIPLY), angularMomentum, EquationElement(.RIGHT), EquationElement(.RIGHT), EquationElement(.MULTIPLY),EquationElement(.LEFT),annihilationOperator(operatorSize), EquationElement(.ADD), creationOperator(operatorSize), EquationElement(.RIGHT)])
}

func positionOperatorQ(operatorSize: Int) -> EquationElement{
    return implimentElementArraySimplification([EquationElement(.SQUAREROOT), EquationElement(.LEFT), EquationElement(.LEFT), mass, EquationElement(.MULTIPLY), angularMomentum, EquationElement(.RIGHT), EquationElement(.DIVIDE), hBar, EquationElement(.RIGHT), EquationElement(.MULTIPLY), positionOperator(operatorSize)])
}

func momentumOperator(operatorSize: Int) -> EquationElement{
    return implimentElementArraySimplification([EquationElement(.SQUAREROOT), EquationElement(.LEFT), EquationElement(.LEFT), mass, EquationElement(.MULTIPLY), angularMomentum, EquationElement(.MULTIPLY), hBar, EquationElement(.RIGHT), EquationElement(.DIVIDE), EquationElement(2.0), EquationElement(.RIGHT), EquationElement(.MULTIPLY), EquationElement(.LEFT), creationOperator(operatorSize), EquationElement(.SUBTRACT), annihilationOperator(operatorSize), EquationElement(.RIGHT)])
}
/*
func momentumOperator(operatorSize: Int) -> EquationElement{
    return implimentElementArraySimplification([i, EquationElement(.MULTIPLY), EquationElement(.SQUAREROOT), EquationElement(.LEFT), EquationElement(.LEFT), mass, EquationElement(.MULTIPLY), angularMomentum, EquationElement(.MULTIPLY), hBar, EquationElement(.RIGHT), EquationElement(.DIVIDE), EquationElement(2.0), EquationElement(.RIGHT), EquationElement(.MULTIPLY), EquationElement(.LEFT), creationOperator(operatorSize), EquationElement(.SUBTRACT), annihilationOperator(operatorSize), EquationElement(.RIGHT)])
}
*/

func momentumOperatorP(operatorSize: Int) -> EquationElement{
    return implimentElementArraySimplification([EquationElement(.SQUAREROOT), EquationElement(.LEFT), EquationElement(1.0), EquationElement(.DIVIDE), EquationElement(.LEFT), hBar, EquationElement(.MULTIPLY), mass, EquationElement(.MULTIPLY), angularMomentum, EquationElement(.RIGHT), EquationElement(.RIGHT), EquationElement(.MULTIPLY), momentumOperator(operatorSize)])
}


func plotableEigenValues(elementArray: [EquationElement], startValue: Float, incrementValue: Float, numberOfSteps: Float, variableIndex: Int) -> [[Float]]{
    var eigenValuesToPlot: [[Float]] = []
    
    for(var loopNumber: Float = 0; loopNumber < numberOfSteps; loopNumber += 1){
        let variableElement: EquationElement = EquationElement(startValue + (loopNumber * incrementValue))
        var finalElementArray = elementArray
        finalElementArray.insert(variableElement, atIndex: variableIndex)
        let element = implimentElementArraySimplification(finalElementArray)
        //print(element.matrix)
        //print("\n\n")
        var eigenValueArray = eigenValues(element)
        sort(&eigenValueArray)
        
        if(loopNumber == 0){
            for(var index = 0; index < eigenValueArray.count; index += 1){
                eigenValuesToPlot.append([])
                //print(index)
            }
        }
        
        for(var index = 0; index < eigenValueArray.count; index += 1){
            eigenValuesToPlot[index].append(eigenValueArray[index])
        }
    }
    //print(eigenValuesToPlot)
    //print("\n\n")
    return eigenValuesToPlot
}

func findMinValueOfPair(a: Float, b: Float) -> Float{
    var result = a
    if(b < a){
        result = b
    }
    return result
}

func cosFAPlusADagger(operatorSize: Int, f: Float) -> EquationElement{
    var result = EquationElement()
    result.type = .MATRIX
    for(var N: Float = 0; N < Float(operatorSize); N += 1){
        for(var M: Float = 0; M < Float(operatorSize); M += 1){
            let partA = implimentElementArraySimplification([EquationElement(.LEFT), EquationElement(N), EquationElement(.FACTORIAL), EquationElement(.MULTIPLY),EquationElement(M), EquationElement(.FACTORIAL), EquationElement(.RIGHT), EquationElement(.POWER), EquationElement(0.5), EquationElement(.MULTIPLY), EquationElement(.EXP), EquationElement(.LEFT),EquationElement(.LEFT), EquationElement(.NEGATE), EquationElement(f), EquationElement(.SQUARED), EquationElement(.RIGHT), EquationElement(.DIVIDE), EquationElement(2), EquationElement(.RIGHT)])
            
            var partB = EquationElement()
            for(var n: Float = 0; n < findMinValueOfPair(N, M); n += 1){
                let partB2 = implimentElementArraySimplification([EquationElement(.LEFT), EquationElement(f), EquationElement(.POWER), EquationElement(.LEFT), EquationElement(N), EquationElement(.ADD), EquationElement(M), EquationElement(.SUBTRACT), EquationElement(2.0), EquationElement(.MULTIPLY), EquationElement(n), EquationElement(.RIGHT), EquationElement(.MULTIPLY), EquationElement(-1), EquationElement(.POWER), EquationElement(.LEFT), EquationElement(.LEFT), EquationElement(N), EquationElement(.ADD), EquationElement(M), EquationElement(.SUBTRACT), EquationElement(2.0), EquationElement(.MULTIPLY), EquationElement(n), EquationElement(.RIGHT), EquationElement(.DIVIDE), EquationElement(2.0), EquationElement(.RIGHT), EquationElement(.MULTIPLY), EquationElement(.LEFT), EquationElement(1.0), EquationElement(.ADD), EquationElement(-1), EquationElement(.POWER), EquationElement(.LEFT), EquationElement(N), EquationElement(.ADD), EquationElement(M), EquationElement(.SUBTRACT), EquationElement(2.0), EquationElement(.MULTIPLY), EquationElement(n), EquationElement(.RIGHT), EquationElement(.RIGHT), EquationElement(.RIGHT), EquationElement(.DIVIDE), EquationElement(.LEFT), EquationElement(2.0), EquationElement(.MULTIPLY), EquationElement(n), EquationElement(.FACTORIAL), EquationElement(.MULTIPLY), EquationElement(.LEFT), EquationElement(N), EquationElement(.SUBTRACT), EquationElement(n), EquationElement(.RIGHT), EquationElement(.FACTORIAL), EquationElement(.MULTIPLY), EquationElement(.LEFT), EquationElement(M), EquationElement(.SUBTRACT), EquationElement(n), EquationElement(.RIGHT), EquationElement(.FACTORIAL), EquationElement(.RIGHT)])
                if(partB.type == .NOTSET){
                    partB = partB2
                }
                else{
                    partB = implimentElementArraySimplification([partB, EquationElement(.ADD), partB2])
                }
            }
            result.matrix.append(implimentElementArraySimplification([partA, EquationElement(.MULTIPLY), partB]).constant)
        }
    }
    //print(result.matrix)
    return result
}

func sinFAPlusADagger(operatorSize: Int, f: Float) -> EquationElement{
    var result = EquationElement()
    result.type = .MATRIX
    for(var N: Float = 0; N < Float(operatorSize); N += 1){
        for(var M: Float = 0; M < Float(operatorSize); M += 1){
            let partA = implimentElementArraySimplification([EquationElement(.LEFT), EquationElement(N), EquationElement(.FACTORIAL), EquationElement(.MULTIPLY),EquationElement(M), EquationElement(.FACTORIAL), EquationElement(.RIGHT), EquationElement(.POWER), EquationElement(0.5), EquationElement(.MULTIPLY), EquationElement(.EXP), EquationElement(.LEFT),EquationElement(.LEFT), EquationElement(.NEGATE), EquationElement(f), EquationElement(.SQUARED), EquationElement(.RIGHT), EquationElement(.DIVIDE), EquationElement(2), EquationElement(.RIGHT)])
            
            var partB = EquationElement()
            for(var n: Float = 0; n < findMinValueOfPair(N, M); n += 1){
                let partB2 = implimentElementArraySimplification([EquationElement(.LEFT), EquationElement(f), EquationElement(.POWER), EquationElement(.LEFT), EquationElement(N), EquationElement(.ADD), EquationElement(M), EquationElement(.SUBTRACT), EquationElement(2.0), EquationElement(.MULTIPLY), EquationElement(n), EquationElement(.RIGHT), EquationElement(.MULTIPLY), EquationElement(-1), EquationElement(.POWER), EquationElement(.LEFT), EquationElement(.LEFT), EquationElement(N), EquationElement(.ADD), EquationElement(M), EquationElement(.SUBTRACT), EquationElement(2.0), EquationElement(.MULTIPLY), EquationElement(n), EquationElement(.SUBTRACT), EquationElement(1.0), EquationElement(.RIGHT), EquationElement(.DIVIDE), EquationElement(2.0), EquationElement(.RIGHT), EquationElement(.MULTIPLY), EquationElement(.LEFT), EquationElement(1.0), EquationElement(.SUBTRACT), EquationElement(-1), EquationElement(.POWER), EquationElement(.LEFT), EquationElement(N), EquationElement(.ADD), EquationElement(M), EquationElement(.SUBTRACT), EquationElement(2.0), EquationElement(.MULTIPLY), EquationElement(n), EquationElement(.RIGHT), EquationElement(.RIGHT), EquationElement(.RIGHT), EquationElement(.DIVIDE), EquationElement(.LEFT), EquationElement(2.0), EquationElement(.MULTIPLY), EquationElement(n), EquationElement(.FACTORIAL), EquationElement(.MULTIPLY), EquationElement(.LEFT), EquationElement(N), EquationElement(.SUBTRACT), EquationElement(n), EquationElement(.RIGHT), EquationElement(.FACTORIAL), EquationElement(.MULTIPLY), EquationElement(.LEFT), EquationElement(M), EquationElement(.SUBTRACT), EquationElement(n), EquationElement(.RIGHT), EquationElement(.FACTORIAL), EquationElement(.RIGHT)])
                if(partB.type == .NOTSET){
                    partB = partB2
                }
                else{
                    partB = implimentElementArraySimplification([partB, EquationElement(.ADD), partB2])
                }
            }
            result.matrix.append(implimentElementArraySimplification([partA, EquationElement(.MULTIPLY), partB]).constant)
        }
    }
    return result
}


func cosq(operatorSize: Int, alpha: Float, beta: Float) -> EquationElement{
    /*
    let f = implimentElementArraySimplification([EquationElement(alpha), EquationElement(.MULTIPLY), EquationElement(.SQUAREROOT), EquationElement(.LEFT), EquationElement(.LEFT), EquationElement(2.0), EquationElement(.MULTIPLY), mass, EquationElement(.MULTIPLY), angularMomentum, EquationElement(.RIGHT), EquationElement(.DIVIDE), hBar, EquationElement(.RIGHT)]).constant
    */
    let f = alpha
    let g = alpha * beta
    
    return implimentElementArraySimplification([cosFAPlusADagger(operatorSize, f), EquationElement(.MULTIPLY), EquationElement(cos(g)), EquationElement(.SUBTRACT), sinFAPlusADagger(operatorSize, f), EquationElement(.MULTIPLY), EquationElement(sin(g))])
}





