参考自ChemAxon的介绍。英文版在[这里](https://docs.chemaxon.com/display/docs/Extended+Connectivity+Fingerprint+ECFP)。

细节部分与rdkit略有不同，但原理基本一致。

# 一、简介

扩展连通性指纹（Extended-Connectivity Fingerprints，ECFPs）是一种圆形拓扑指纹，可用于分子表示、相似性搜索、构效关系建模。它在药物设计中有诸多应用，是最受欢迎的相似性比对工具。

# 二、应用

ECFPs是非常高效且流行的搜索工具，并有着广泛的应用。

ECFPs最早应用在高通量筛选（high-throughput screening，HTS）中，用于分析HTS中假阳性和假阴性先导物。此外，ECFPs也经常应用在基于配体的虚拟筛选中，用来判断化合物是否有活性。大量的实验结果表明，这种圆形拓扑指纹是效果最好的搜索工具之一。

在药物研究中，许多领域都涉及到相似性搜索，都用到了分子结构中蕴含的大量信息，例如化合物聚类，以及化合物库分析等。

除了相似性搜索，ECFPs也可以用来识别是否含有一个特定的子结构。因此，它也经常应用在QSAR和QSPR建模中（例如ADMET属性预测），并优化先导物性质。

# 三、特性

ECFPs的主要属性有如下几点：
* 它通过圆形原子邻域（circular atom neighborhoods）来表示分子结构
* 计算速度快
* 这些特征代表特定子结构的存在与否
* 没有经过预定义，可以表示大量不同的分子特性（包括立体化血信息）
* 它用来表示功能基团（functionality）的是否存在，对于分析分子活性至关重要
* 它的生成算法非常灵活，可以针对不同的场景，生成不同的指纹

# 四、与路径指纹的对比

基于路径的指纹（path-based fingerprints）广泛用于子结构搜索中预过滤（pre-filtering）。与这种指纹相反，ECFPs并不适合子结构搜索，而适合快速高效地筛选整个分子和进行相似性搜索。ECFPs可以为相似性搜索提供充分的结果，更接近药物化学家的需求。

# 五、表示方法

ECFPs在ChemAxon中有如下两种表示方式。

## 1.整数标识列（list of integer identifiers）

* ECFPs可以通过一组不定长的整数标识符来表示，这是最原始和最精确的表示方法。每个标识符代表了一个特定的子结构，更准确地说是分子里的一种圆形原子邻域。这种整数标识符是按升序排列。
* 这些标识符也可以被理解为是一个超大的虚拟比特串的索引。比特串中的每个位置表示一个特定子结构的特征。因为这个虚拟的比特串非常庞大且稀疏，所以它并没有被显式地存储，而是用代表比特串的索引来组成这个不定长的整数列表。由于技术原因，这种特征标识符存储为一种有符号的数值，因此是可正可负的。
* 默认情况下，在一个整数标识列中，每种标识符最多仅含一个（rdkit默认不会去重）。然而，在一些特定的场景下，也需要考虑某种标识符在分子中出现的频数，这种ECFPs的变体叫做ECFC。在ChemAxon的工具中，有一个参数控制是否统计频数，默认模式为不统计，统计时称为ECFC模式。

## 2.定长比特串（fixed-length bit string）

* 传统的二进制分子表示方法是通过固定长度的比特串实现的。而ECFPs也可以通过“折叠”操作来形成定长比特串，也就是将底层的整数标识列压缩成一个更短长度的定长比特串（就是一个向量，常见的长度有1024,2048等）。
* 与上面的整数标识列表相比，ECFPs的定长比特串简化了比对和相似性计算，也相应减少了存储开销，尤其是大分子。另一方面，压缩的操作可能会引入冲突，因为压缩后，多个不同的子结构特征可能会用同一个比特位表示。这种情况下，会丢失一些分子信息，造成特征质量和可解释性的下降。
* 需要注意的是，由标识列可以得到定长比特串，但没法从定长比特串得到标识列。定长比特串可以看作是标识列的有损压缩。

# 六、生成流程

可基于ChemAxon的GenerateMD命令行，或API生成，感兴趣的可以自行查看。

## 1.初始化原子标识符

* ECFP生成时，首先对给定分子的每一个非氢原子分配一个初始整数标识符。该标识符是通过把原子的属性（例如原子序号、连接数量等）打包，再经哈希函数转变成一个整数变量而得到，包含了相应原子的化学信息。
* 需要考虑的原子属性是ECFPs中一个重要的参数，可以完全由自己设置。

## 2.标识符的迭代更新

* 初始化后，会进行一系列迭代操作，将某一原子的初始标识符与邻近原子的标识符合并，直到到达设置的半径为止。每一轮迭代都会捕捉距中心原子越来越远的原子信息，最终经哈希运算，编码成为一个整数值，这些整数值合并形成一个整数标识列。

<center>

![1](https://github.com/dreadlesss/extended_connectivity_fingerprints/blob/master/png/1.png)
Fig. 1. 对某一特定原子进行迭代更新的示意图</center>

* 这种迭代更新的过程基于一种摩根算法（Morgan algorithm）。

## 3.标识符去重

* 在最终生成指纹的过程中，会去除等价的原子标识符。等价意味着他们包含一组相同的键，或他们有相同的哈希值。在这步操作中，如果标识符的频数需要保留（ECFC模式），相同标识符出现的次数将会被保存。
* 下面两张图片展示了生成一个分子ECFPs的完整流程

<center>

![2](https://github.com/dreadlesss/extended_connectivity_fingerprints/blob/master/png/2.png)
Fig. 2. 从分子到整数标识列的生成过程

![3](https://github.com/dreadlesss/extended_connectivity_fingerprints/blob/master/png/3.png)
Fig. 3. 从整数标识列到定长比特串的生成过程</center>

# 七、参数设置

ChemAxon中的ECFPs有一系列的参数，可以通过XML来进行设置。

## 1.主要参数

主要参数有3个：ECFPs的最大直径、指纹长度、标识符频数。

* 直径

这个参数规定了每个原子要考虑的最大圆形邻域。默认直径是4。这是个很重要的参数，会影响子结构的数量和大小，从而影响整数标识列的长度，也会影响到定长比特串的表示。

ECFPs可以通过该参数进行区分。例如，ECFP_4表示直径为4，ECFP_6表示直径为6（rdkit中的参数为半径radius，当radius=2时，与ECFP_4等价）。

根据Rogers and Hahn的报道，直径如何取值取决于具体需求。通常直径为4时，可以满足相似性搜索和聚类需求。对于活性预测，有大结构片段通常会比较有利，因此可以尝试6或8。

* 长度

这个参数指定了定长比特串的长度。默认是1024位。

当增加长度时，会减少冲突的几率，从而减少信息的损失，但需要花更多的空间和时间来存储和计算。

* 频数

这个参数控制是否要记录每个整数标识符出现的总次数。默认是不记录，即每个子结构只存储1次。


## 2.原子属性

在ECFPs中，原子属性集也是一个重要的参数，决定了原子的特性。不同的属性会产生不同的圆形指纹，可用于不同的需求中。

默认的ECFP参数适用于大多数的应用场景，尤其在基于子结构特异信息的相似性搜索中。默认情况下，对每个原子的以下信息会被考虑进来：

* 原子序号
* 中心原子的相邻重原子（非氢原子）数量
* 中心原子的相邻氢原子数量
* 形式电荷
* 附加属性：该原子是否是环中的一部分

其他的内置属性还有：质量信息、连接信息、键信息、立体信息等。如果想修改这些配置，可以在XML中自己定义。

# 八、功能基指纹

ECFP默认的标识符包含了高度特异的原子信息，可以表示大量的子结构特征。然而，在不同的场景中，可能需要不同的抽象方式。例如，在环上某一位置进行氯取代或溴取代时，它们的功能基本是等价的，但是ECFP认为它们是两种不同的结构。另一个使用圆形指纹的新方向就是归纳每个原子的药效团，也就是将ECFPs转换为拓扑药效团指纹。

这些变体，叫做功能基指纹（Functional-Class Fingerprints，FCFPS），它们对ECFP进行了泛化处理，更关注原子在功能上所发挥的作用，而不关注是否是特异的。

ChemAxon也支持FCFPs，可通过XML来自行定义。

# 九、特征可视化

JChem提供了查看ECFP指纹结构的功能。ECFP可以表示为整数标识列，或定长比特串，这些特征中的某一位代表了输入分子中一个特定的子结构。ECFPFeatureLookup类可以根据给定的整数标识和某个位来提取出子结构。

在某些应用场景中，需要对特定化合物的指纹进行对比，以便得出有哪些特征很重要。这时，展示出真实的子结构就会非常直观地看到结果。

# 十、rdkit计算ECFPs

在rdkit中指纹使用位向量[Bit Vector](http://www.rdkit.org/docs/GettingStartedInPython.html#bit-vectors)来表示，位向量是一种可以高效存放二进制值的容器。rdkit中包含了两种内部存储方式不同的指纹，分别是SparseBitVects和ExplicitBitVects。这两种类型可以轻易地进行转换，用于不同的目的。

* SparseBitVects，是一种非常大，非常稀疏的向量，在该类型上进行索引和提取会非常快，有点类似于list of integer identifiers。
* ExplicitBitVects，记录某一位上是否存在，计算速度通常比SparseBitVects快，但会占用更多内存，类似于fixed-length bit string。
* 生成SparseBitVects，需要设置参数mol和半径radius
```python
>>> from rdkit.Chem import AllChem
>>> mol = Chem.MolFromSmiles('Cc1ccccc1')
>>> AllChem.GetMorganFingerprint(mol, radius=2)
<rdkit.DataStructs.cDataStructs.UIntSparseIntVect at 0x7fbbe9654ee0>
```
* 生成ExplicitBitVects，除了mol和radius，还要传入向量长度nBits
```python
>>> AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)
<rdkit.DataStructs.cDataStructs.ExplicitBitVect at 0x7fbbe9493490>
```