-- MariaDB dump 10.19  Distrib 10.5.22-MariaDB, for Linux (x86_64)
--
-- Host: localhost    Database: hg38
-- ------------------------------------------------------
-- Server version	10.5.22-MariaDB

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8mb4 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `ncbiRefSeqLink`
--

DROP TABLE IF EXISTS `ncbiRefSeqLink`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `ncbiRefSeqLink` (
  `id` varchar(255) NOT NULL,
  `status` varchar(255) NOT NULL,
  `name` varchar(255) NOT NULL,
  `product` varchar(255) NOT NULL,
  `mrnaAcc` varchar(255) NOT NULL,
  `protAcc` varchar(255) NOT NULL,
  `locusLinkId` varchar(255) NOT NULL,
  `omimId` varchar(255) NOT NULL,
  `hgnc` varchar(255) NOT NULL,
  `genbank` varchar(255) NOT NULL,
  `pseudo` varchar(255) NOT NULL,
  `gbkey` varchar(255) NOT NULL,
  `source` varchar(255) NOT NULL,
  `gene_biotype` varchar(255) NOT NULL,
  `gene_synonym` varchar(255) NOT NULL,
  `ncrna_class` varchar(255) NOT NULL,
  `note` longblob NOT NULL,
  `description` longblob NOT NULL,
  `externalId` varchar(255) NOT NULL,
  PRIMARY KEY (`id`),
  KEY `name` (`name`),
  KEY `mrnaAcc` (`mrnaAcc`),
  KEY `protAcc` (`protAcc`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_swedish_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2025-08-17  4:31:57
